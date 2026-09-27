#!/usr/bin/env python3
"""Replace classes of a panel with the same classes from builds that masked only their own copies.

    merge_own.py MAIN.tsv.gz OUT.tsv.gz CLASS=PANEL.tsv.gz ...

Every build was made from the same manifest (so the class ids agree and k-mers shared between two
units are gone from both). In MAIN the background masks are the union of all classes' loci, so a
k-mer that one class shares with a paralog copy of ANOTHER class's locus (MW2/MW3 inside the opsin
array for OPN1LW, the CDY2 loci for CDY1) is never background. In CLASS's own build only its own
copies are masked and such k-mers drop. Output is gzip with a zero timestamp.
"""
import gzip
import io
import sys


def read(path):
    head, classes, rows = [], {}, []
    for line in gzip.open(path, "rt"):
        line = line.rstrip("\n")
        if line.startswith("##class\t"):
            f = dict(x.split("=", 1) for x in line.split("\t")[1:])
            classes[f["id"]] = (f["name"], line)
        elif line.startswith("#"):
            head.append(line)
        elif line:
            rows.append(line)
    return head, classes, rows


def main():
    main_path, out_path, specs = sys.argv[1], sys.argv[2], sys.argv[3:]
    head, classes, rows = read(main_path)
    ids = {name: i for i, (name, _) in classes.items()}
    lines = {i: line for i, (_, line) in classes.items()}
    take = {}
    for spec in specs:
        name, path = spec.split("=", 1)
        h, c, r = read(path)
        i = ids[name]
        assert c[i][0] == name, f"{path}: class id {i} is {c[i][0]}, not {name} (manifests differ?)"
        lines[i] = c[i][1]
        take[i] = [x for x in r if x.split("\t", 2)[1] == i]
    kept = [x for x in rows if x.split("\t", 2)[1] not in take]
    for i in take:
        kept += take[i]
    kept.sort(key=lambda x: (int(x.split("\t")[1]), int(x.split("\t")[2])))
    names = ", ".join(classes[i][0] for i in take)
    buf = io.StringIO()
    for h in head:
        if h.startswith("#kmer"):
            buf.write(f"##filter=merge_own.py: {names} from builds whose masks hold only the class's own copies "
                      "(masks/*.own_copies.bed)\n")
            for i in sorted(lines, key=int):
                buf.write(lines[i] + "\n")
        buf.write(h + "\n")
    for x in kept:
        buf.write(x + "\n")
    with open(out_path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as z:
        z.write(buf.getvalue().encode())
    for i in sorted(take, key=int):
        print(f"[merge_own] {classes[i][0]}: {len(take[i])} k-mers from its own-copy build", file=sys.stderr)


if __name__ == "__main__":
    main()
