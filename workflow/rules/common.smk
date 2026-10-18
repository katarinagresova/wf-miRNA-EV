# Functions of the rule files. They are in a file without rules, as the Snakemake linter recommends.

import re
import subprocess


def git_commit():
    """Commit of the workflow source, or 'unknown' (for example, a release archive)."""
    out = subprocess.run(["git", "-C", str(workflow.current_basedir), "describe", "--always", "--dirty", "--tags"],
                         capture_output=True, text=True)
    return out.stdout.strip() if out.returncode == 0 else "unknown"


def ref_source(source, download):
    """Path of a reference file: the download path if the config gives a URL (or a URL and md5), else the local path
    in the config."""
    url = source["url"] if isinstance(source, dict) else source
    return download if re.match(r"(https?|ftp)://", url) else url


def ncbi_source(acc, ext, download):
    """Path of an NCBI record (rDNA or contaminant): the local copy in refs.ncbi_dir, or the download path."""
    return f"{REFS['ncbi_dir']}/{acc}.{ext}" if REFS["ncbi_dir"] else download
