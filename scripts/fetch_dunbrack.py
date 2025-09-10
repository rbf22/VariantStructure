import tarfile
import pathlib
import subprocess


def main():
    root = pathlib.Path(__file__).resolve().parents[2]
    tarball = root / "data" / "bbdep02.May2010.tar.gz"

    # Ensure Git LFS has pulled the file
    subprocess.run(["git", "lfs", "pull"], check=True)

    # Extract
    outdir = root / "data" / "dunbrack"
    outdir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tarball, "r:gz") as tar:
        tar.extractall(outdir)

    print(f"Dunbrack library extracted to {outdir}")
