import re
import sys
import os


def trim(docstring):
    """
    Trim the docstring (trim function from PEP-257)
    :param docstring: str - docstring for trimming
    :return: str - trimmed docstring
    """
    if not docstring:
        return ""
    # Convert tabs to spaces (following the normal Python rules)
    # and split into a list of lines:
    lines = docstring.expandtabs().splitlines()
    # Determine minimum indentation (first line doesn't count):
    indent = sys.maxsize
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indent = min(indent, len(line) - len(stripped))
    # Remove indentation (first line is special):
    trimmed = [lines[0].strip()]
    if indent < sys.maxsize:
        for line in lines[1:]:
            trimmed.append(line[indent:].rstrip())
    # Strip off trailing and leading blank lines:
    while trimmed and not trimmed[-1]:
        trimmed.pop()
    while trimmed and not trimmed[0]:
        trimmed.pop(0)

    # Return a single string:
    return "\n".join(trimmed)


def parse_docstring_snakemake(docstring):
    """
    Parse snakemake docstring.
    :param docstring: str - contents of a snakemake docstring
    :return: dict - dict with ["description", "params", "inputs", "outputs"] keys filled out
    """

    # prepare defaults:
    description = ""
    params = []
    inputs = []
    outputs = []

    # prepare regexps
    keywords = ['params?', 'input', 'output']
    end_keyword = '(?:(?=:' + ')|(?=:'.join(keywords) + ')|\Z|""")'

    keywords_regex = re.compile(":(?:" + "|".join(keywords) + ")")
    param_regex = re.compile(":param (?P<name>[\w]+): (?P<doc>.*?)" + end_keyword, re.DOTALL)
    input_regex = re.compile(":input (?P<name>[\w]+): (?P<doc>.*?)" + end_keyword, re.DOTALL)
    output_regex = re.compile(":output (?P<name>[\w]+): (?P<doc>.*?)" + end_keyword, re.DOTALL)

    if docstring is not None:

        # trim docstring
        docstring = trim(docstring)
        params_returns_desc = None

        # try to find the first keyword
        match = keywords_regex.search(docstring)
        if match:
            desc_end = match.start()
            params_returns_desc = docstring[desc_end:].strip()
            description = docstring[:desc_end].rstrip()
        else:
            description = docstring.rstrip()

        if params_returns_desc:
            # get params:
            params = [{"name": name, "doc": trim(doc).strip()} for name, doc in param_regex.findall(params_returns_desc)]
            # get inputs:
            inputs = [{"name": name, "doc": trim(doc).strip()} for name, doc in input_regex.findall(params_returns_desc)]
            # get outputs:
            outputs = [{"name": name, "doc": trim(doc).strip()} for name, doc in output_regex.findall(params_returns_desc)]

    return {
        "description": description,
        "params": params,
        "inputs": inputs,
        "outputs": outputs
    }


def parse_snakemake_str(snakemake):
    """
    Parse the snake file and return rule names and their docstrings.
    :param snakemake: str - contents of a Snakefile
    :return: list(str, str), list(str) - list of all rules with docstrings and list of all rules
    """
    # define keywords
    keywords = ['rule .*:', 'input:', 'output:', 'params:', 'threads:', 'run:', 'shell:', 'log:', 'benchmark:', 'message:', 'resources:', 'version:']
    # keywords_regex = re.compile("(?:" + "|".join(keywords) + ")")
    end_keyword = '(?:(?=' + ')|(?='.join(keywords) + ')|\Z)'

    # define regexps for rules and docstrings and for rules only
    rule_docs_regex = re.compile('.*?rule (?P<rulename>[\w]+):.*?"""(?P<docstring>.*?)""".*?' + end_keyword, re.DOTALL)
    rule_all_regex = re.compile('.*?rule (?P<rulename>[\w]+):.*?' + end_keyword, re.DOTALL)

    # return all rules with docstrings and all rules
    return rule_docs_regex.findall(snakemake), rule_all_regex.findall(snakemake)


def parse_snakemake(filename):
    """
    Parse the snake file and return rule names and their docstrings.
    :param filename: str - filename of the Snakefile
    :return: list(str, str), list(str) - list of all rules with docstrings and list of all rules
    """
    if type(filename) is str:
        with open(filename) as f:
            return parse_snakemake_str(f.read())


def generate_rst(snakemake_file, docs_file, append=False):
    """
    Generate .rst doc file from the docstring.
    :param snakemake_file: str - snakemake filename (with rules and docstrings)
    :param docs_file: str - filename, where to output documentation
    :param append: bool - append to docs_file?
    :return: int, int - number of written rules with docstrings, number of rules with missing docstrings
    """
    # prepare variables
    written_rules = 0
    rules_wo_docs = 0

    # templates:
    rule_template = """rule {name}\n{underscore}\nlocated in: '<SnakeLines_dir>/{filename}'\n\n{description}\n\n"""
    param_begin_template = "{pname}:\n"
    param_template = "\t{name}: {doc}\n"
    rule_end_template = "\n"
    undoc_rules_beg_template = "undocumented rules\n------------------\nWARNING: found {many:2d} undocumented rules:\n"
    undoc_rules_template = "\t- rule {rule} is UNDOCUMENTED\n"

    # print("Parsing", snakemake_file, docs_file)

    # try to convert:
    try:
        # parse snakefile:
        docs_rules, all_rules = parse_snakemake(snakemake_file)

        # setup defaults
        names = []
        rst = ""

        # continue only if some rules found
        if len(docs_rules) > 0:

            # undocumented rules:
            names, docstrings = list(zip(*docs_rules))

            # write rst
            for name, docstring in zip(names, docstrings):

                # parse docstring:
                docstr_dict = parse_docstring_snakemake(docstring)

                # write begining
                underscore = '-' * (len(name) + 5)
                rst += rule_template.format(name=name, description=docstr_dict['description'], underscore=underscore, filename=snakemake_file)

                # write all inputs/outputs/params
                rst += param_begin_template.format(pname='Input(s)')
                for param in docstr_dict['inputs']:
                    rst += param_template.format(name=param['name'], doc=param['doc'])
                if len(docstr_dict['outputs']) > 0:
                    rst += param_begin_template.format(pname='Output(s)')
                    for param in docstr_dict['outputs']:
                        rst += param_template.format(name=param['name'], doc=param['doc'])
                if len(docstr_dict['params']) > 0:
                    rst += param_begin_template.format(pname='Param(s)')
                    for param in docstr_dict['params']:
                        rst += param_template.format(name=param['name'], doc=param['doc'])

                # finalize
                rst += rule_end_template.format()
                written_rules += 1

        # write undocumented rule names
        if len(all_rules) - len(docs_rules) > 0:

            # extract undocumented rule names
            undoc_rules = list(filter(lambda x: x not in names, all_rules))
            rules_wo_docs = len(undoc_rules)

            # write undocumented rules
            rst += undoc_rules_beg_template.format(many=rules_wo_docs)
            for undoc_rule in undoc_rules:
                rst += undoc_rules_template.format(rule=undoc_rule)

        # write the rst into a file:
        if rst != "":
            # open and write file
            dirname = os.path.dirname(docs_file)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            with open(docs_file, 'a' if append else 'w') as f:
                f.write(rst)

    except Exception as e:

        # write what is wrong
        print(repr(e))

        # set defaults
        written_rules = 0
        rules_wo_docs = 0

    return written_rules, rules_wo_docs


def crawl_snakedir(input_snakedir, docs_snakedir, depth, first_line=""):
    """
    Crawl input snakedir, if found Snakefile, create docs and output them into docs dir.
    :param input_snakedir: str - path to snakedir to crawl
    :param docs_snakedir: str - path to snakedir, where to output documentation
    :param depth: int - depth of recursion
    :param first_line: str - first line of the .rst string
    :return: int, int, int, str - number of files (re)created, number of written rules with docstrings, number of rules with missing docstrings, rst string generated so far
    """
    files_done = 0
    rules_written = 0
    rules_wo_docs = 0
    summary_name = 'summary.rst'

    if not input_snakedir.endswith('/'):
        input_snakedir += '/'

    if not docs_snakedir.endswith('/') and docs_snakedir != "":
        docs_snakedir += '/'

    # define templates
    rst = ""
    dir_template = "{depth}- {filename} ({doc_rules}/{all_rules} documented/all rule(s), {file_count} file(s))\n"
    first_line_template = "{first_line}\n{underscore}\nThis document contains {doc_rules}/{all_rules} documented/all rule(s) in {file_count} file(s):\n"

    # look at the dir:
    root, dirs, files = next(os.walk(input_snakedir))
    files_filtered = list(filter(lambda filename: filename.endswith('.snake') or filename.startswith('Snakefile'), files))

    # first do files
    rst_path = docs_snakedir + summary_name
    relative_rst_path = '/'.join(rst_path.split('/')[2:])
    for filename in files_filtered:
        rw, rwo_docs = generate_rst(input_snakedir + filename, rst_path, append=files_done > 0)
        if rw > 0 or rwo_docs > 0:
            files_done += 1
            rules_written += rw
            rules_wo_docs += rwo_docs

    if len(files_filtered) > 0:
        rst += '{depth}.. toctree::\n'.format(depth='\t'*depth)
        rst += '{depth}{filename}\n\n'.format(depth='\t'*(depth+1), filename=relative_rst_path)

    # recursive do directories
    for directory in dirs:
        fd, rw, rwo_docs, rst_so_far = crawl_snakedir(input_snakedir + directory, docs_snakedir + directory, depth + 1)
        files_done += fd
        rules_written += rw
        rules_wo_docs += rwo_docs

        if fd > 0:
            rst += dir_template.format(filename=directory, filepath=docs_snakedir + directory + '/' + summary_name,
                                       doc_rules=rw, all_rules=rw + rules_wo_docs, file_count=fd, depth="\t" * depth)
            rst += rst_so_far

    # add first line:
    if depth == 1 and first_line != "":
        underscore = "=" * len(first_line)
        rst = first_line_template.format(first_line=first_line, underscore=underscore, doc_rules=rules_written, all_rules=rules_written + rules_wo_docs, file_count=files_done) + rst + '\n'

    return files_done, rules_written, rules_wo_docs, rst


def write_rst(rst, docs_snakedir):
    """
    Write the outline rst to file.
    :param rst: str - generated documentation in string
    :param docs_snakedir: str - directory of the documentation
    :return: None
    """
    # setup
    summary_name = 'outline.rst'

    if not docs_snakedir.endswith('/'):
        docs_snakedir += '/'

    # write the summary doc:
    if rst != "":
        # open and write file
        dirname = os.path.dirname(docs_snakedir + summary_name)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(docs_snakedir + summary_name, 'w') as f:
            f.write(rst)


if __name__ == '__main__':
    # run documentation generation for rules and pipelines:
    files, rules, missing, outline = crawl_snakedir('rules',  'docs/rules', 1, "Rules")
    files_p, rules_p, missing_p, outline_p = crawl_snakedir('pipeline', 'docs/pipeline', 1, "Pipeline")

    # write the summary doc:
    write_rst(outline, 'docs/rules')
    write_rst(outline_p, 'docs/pipeline')

    # if we have all pipeline rules documented, then commit else fail
    exit(missing_p)
