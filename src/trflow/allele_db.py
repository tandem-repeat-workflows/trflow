import gzip


def get_info_tags(line):
    sl = line.split()
    format, values = sl[8], sl[9]

    tags = format.split(":")
    vals = values.split(":")

    return {tag: val for tag, val in zip(tags, vals)}


def get_alleles(vcf_path):
    with gzip.open(vcf_path, "rb") as file:
        for line in file:
            line = line.decode("utf8")
            if line[0] == "#":
                continue

            sl = line.split()
            ref, alts = sl[3], sl[4]
            locus = sl[7].split(";")[0].replace("TRID=", "")
            if len(sl) <= 9:
                assert False, sl
                continue
            gt = sl[9].split(":")[0]

            if gt == "0/0":
                assert alts == ".", alts
                alleles = [ref, ref]
            elif gt == "1/2":
                assert "," in alts
                alleles = alts.split(",")
            elif gt == "1/1":
                assert "," not in alts
                alleles = [alts, alts]
            elif gt == "0/1":
                assert "," not in alts
                alleles = [ref, alts]
            elif gt == "1":
                assert "," not in alts
                alleles = [alts]
            elif gt == "0":
                assert alts == "."
                alleles = [ref]
            elif gt == ".":
                continue
            else:
                assert False, gt

            info_tags = get_info_tags(line)

            yield locus, alleles, info_tags


def create_allele_db(manifest, db_path):
    r"""Generate an allele database

    Parameters
    ----------
    manifest : dict
        A dictionary where keys are sample names and values are paths to TRGT VCFs
    db_path : str
        A path to the allele db file

    Examples
    --------
    >>> manifest = {"SampleA": "/path/to/a.vcf.gz", "SampleB": "/path/to/b.vcf.gz"}
    >>> crate_allele_db(manifest, "output/allele_db.gz")
    """

    with gzip.open(db_path, "wb") as file:
        for sample, vcf_path in manifest.items():
            for trid, alleles, tags in get_alleles(vcf_path):
                alleles = ",".join(alleles)
                alid = tags["AL"]
                file.write(f"{trid}\t{sample}\t{alid}\t{alleles}\n".encode("utf-8"))
