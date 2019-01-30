# TODO this should be changed to the version of SnakeLines, e.g. 1.0.2


def get_git_revision_hash(script_dir):
    """
    Get git hash.
    :param script_dir: str - path to script_dir
    :return: (str, str)- hash of git repo and its shorter version
    """
    try:
        print('%s/.git/refs/remotes/origin/master' % script_dir)
        with open('%s/.git/refs/remotes/origin/master' % script_dir) as f:
            print('%s/.git/refs/remotes/origin/master' % script_dir)
            res = f.read().strip()
            short = res[:7]
        return res, short
    except Exception:
        # User does not have snakemake downloaded from git
        return None, None


def store_revision(script_dir, report_dir):
    """
    Get git hash.
    :param script_dir: str - path to script_dir
    :param report_dir: str - output directory for reported files
    """
    long_hash, short_hash = get_git_revision_hash(script_dir)
    if long_hash and short_hash:
        with open('{}/version.txt'.format(report_dir), 'w') as out:
            out.write('Long hash: {}\n'.format(long_hash))
            out.write('Short hash: {}\n'.format(short_hash))