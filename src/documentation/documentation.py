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


def generate_rst(snakemake_file, docs_file):
    """
    Generate .rst doc file from the docstring.
    :param snakemake_file: str - snakemake filename (with rules and docstrings)
    :param docs_file: str - filename, where to output documentation
    :return: int, int -number of written rules with docstrings, number of rules with missing docstrings
    """
    # prepare variables
    written_rules = 0
    rules_wo_docs = 0

    # templates:
    rule_template = """rule {name}\n{underscore}\n{description}\n\n"""
    param_begin_template = "{pname}:\n"
    param_template = "\t{name}: {doc}\n"
    rule_end_template = "\n"
    undoc_rules_beg_template = "WARNING: found {many:2d} undocumented rules:\n"
    undoc_rules_template = "\trule {rule} is UNDOCUMENTED...\n"

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
                rst += rule_template.format(name=name, description=docstr_dict['description'], underscore=underscore)

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
            with open(docs_file, 'w') as f:
                f.write(rst)

    except Exception as e:

        # write what is wrong
        print(repr(e))

        # set defaults
        written_rules = 0
        rules_wo_docs = 0

    return written_rules, rules_wo_docs


def crawl_snakedir(input_snakedir, docs_snakedir):
    """
    Crawl input snakedir, if found Snakefile, create docs and output them into docs dir.
    :param input_snakedir: str - path to snakedir to crawl
    :param docs_snakedir: str - path to snakedir, where to output documentation
    :return: int, int, int - number of files (re)created, number of written rules with docstrings, number of rules with missing docstrings
    """
    files_done = 0
    rules_written = 0
    rules_wo_docs = 0
    summary_name = 'summary.rst'

    if not input_snakedir.endswith('/'):
        input_snakedir += '/'

    if not docs_snakedir.endswith('/'):
        docs_snakedir += '/'

    # define templates
    rst = ""
    file_template = ":download:`File: {filename} <{filepath}>` ({doc_rules}/{all_rules} documented/all rule(s))\n"
    dir_template = ":download:`Directory: {filename} <{filepath}>` ({doc_rules}/{all_rules} documented/all rule(s), {file_count} file(s))\n"

    # look at the dir:
    root, dirs, files = next(os.walk(input_snakedir))

    # first do files
    for filename in files:
        if filename.endswith('.snake') or filename.startswith('Snakefile'):
            rst_path = docs_snakedir + filename + '.rst'
            rst_file = filename + '.rst'
            rw, rwo_docs = generate_rst(input_snakedir + filename, rst_path)
            if rw > 0 or rwo_docs > 0:
                files_done += 1
                rules_written += rw
                rules_wo_docs += rwo_docs
            rst += file_template.format(filename=rst_file, filepath=rst_path, doc_rules=rw, all_rules=rw + rwo_docs)

    # recursive do directories
    for directory in dirs:
        fd, rw, rwo_docs = crawl_snakedir(input_snakedir + directory, docs_snakedir + directory)
        files_done += fd
        rules_written += rw
        rules_wo_docs += rwo_docs

        rst += dir_template.format(filename=directory, filepath=docs_snakedir + directory + '/' + summary_name,
                                   doc_rules=rw, all_rules=rw + rules_wo_docs, file_count=fd)

    # write the summary doc:
    if rst != "":
        # open and write file
        dirname = os.path.dirname(docs_snakedir + summary_name)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        with open(docs_snakedir + summary_name, 'w') as f:
            f.write(rst)

    return files_done, rules_written, rules_wo_docs


if __name__ == '__main__':
    # run documentation generation for rules and pipelines:
    files, rules, missing = crawl_snakedir('rules', 'docs/rules')
    files_p, rules_p, missing_p = crawl_snakedir('pipeline', 'docs/pipeline')

    # if we have all pipeline rules documented, then commit else fail
    exit(missing_p)
