import os
import subprocess
import shutil
import sys

binary = " VeryFastTree"

binaries = {
    "OpenMP": [],
    "OpenMP-AVX2": ["-DUSE_AVX2=ON"],
    "OpenMP-AVX512": ["-DUSE_AVX512=ON"],
}


def main(version):
    for name, options in binaries.items():
        shutil.rmtree("win-build", ignore_errors=True)
        shutil.copytree(os.getcwd(), "win-build")
        os.chdir("win-build")
        subprocess.call(["cmake", "."] + options)
        subprocess.call(["cmake", "--build", "."])
        os.chdir("..")
        os.mkdir("win-build/release")
        for file in os.listdir("win-build/Debug/"):
            if file.endswith(".exe"):
                shutil.copy("win-build/Debug/" + file, "win-build/release/")
        shutil.copy("README.md", "win-build/release")
        shutil.copy("LICENSE", "win-build/release")
        shutil.make_archive(binary + version + "Win64-" + name, 'zip', "win-build/release")
        shutil.rmtree("win-build", ignore_errors=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("version required")
        exit(-1)
    main(sys.argv[1])
