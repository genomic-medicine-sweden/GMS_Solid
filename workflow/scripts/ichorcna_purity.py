import csv

COLUMNS = ["Sample", "Tumor Fraction", "Ploidy", "Subclone Fraction", "GC-Map correction MAD"]


def main(params_txt, output):
    with open(params_txt) as f:
        lines = [line.rstrip("\n") for line in f]

    selected = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    key_values = {}
    for line in lines[2:]:
        if line.startswith("init\t"):
            break
        if ":\t" in line:
            key, value = line.split(":\t", 1)
            key_values.setdefault(key, value)

    row = {
        "Sample": selected["Sample"],
        "Tumor Fraction": selected["Tumor Fraction"],
        "Ploidy": selected["Ploidy"],
        "Subclone Fraction": key_values.get("Subclone Fraction", "NA"),
        "GC-Map correction MAD": key_values.get("GC-Map correction MAD", "NA"),
    }

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main(snakemake.input.params_txt, snakemake.output.purity)
