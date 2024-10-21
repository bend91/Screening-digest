import struct
# import html2text


# def parse(val):

def note_parser(notes):
    parsed_note = {}
    note_offset = 0
    while note_offset < len(notes):
        note_sub = notes[note_offset:]
        header = note_sub[note_sub.find("<")+1:note_sub.find(">")]
        content = note_sub[note_sub.find(">") + 1:note_sub.find("<")]
        parsed_note[header] = content
        note_offset = note_offset + note_sub.find("")
    return parsed_note
    # header1 =
    # content1 = notes[notes[note_offset:].find(">") + 1:notes[note_offset:].find("<")]
    # parsed_note[header1] = content1


def parse_dict(obj):
    if isinstance(obj, dict):
        for key in obj:
            if isinstance(obj[key], dict):
                parse_dict(obj[key])
    return obj


def unpack(fileobject, size, mode):
    return struct.unpack(">" + mode, fileobject.read(size))[0]


def parse_snapgene_file(filepath):
    fileobject = open(filepath, "rb")
    if fileobject.read(1) != b"\t":
        print("Not a snapgene file")
    length = unpack(fileobject, 4, "I")
    title = fileobject.read(8).decode("ascii")
    if title != "SnapGene":
        print("Not a snapgene file")
    data = {
        "isDNA": unpack(fileobject, 2, "H"),
        "exportVersion": unpack(fileobject, 2, "H"),
        "importVersion": unpack(fileobject, 2, "H")
    }
    while True:
        next_byte = fileobject.read(1)
        if next_byte == b"":
            print("End of file")
            break
        block_size = unpack(fileobject, 4, "I")
        if ord(next_byte) == 0:
            props = unpack(fileobject, 1, "b")
            data["dna"] = {
                "topology": "circular" if props & 0x01 else "linear",
                "strandedness": "double" if props & 0x02 > 0 else "single",
                "damMethylated": props & 0x04 > 0,
                "dcmMethylated": props & 0x08 > 0,
                "ecoKIMethylated": props & 0x010 > 0,
                "length": block_size - 1
            }
            data["seq"] = fileobject.read(block_size - 1).decode("ascii")

        elif ord(next_byte) == 6:
            block_content = fileobject.read(block_size).decode("utf-8")
            note_data = parse_dict(block_content)
            # data["notes"] = note_data.get("Notes")
            data["notes"] = note_data

        else:
            fileobject.read(block_size)
            pass
    fileobject.close()
    return data


def main():
    filepath = "/Users/benjamindraper/Library/CloudStorage/Dropbox/01 UCL/Archive/Plasmids/[Addgene]8xNFAT-ZsG-hCD8.dna"
    file_data = parse_snapgene_file(filepath)


if __name__ == "__main__":
    main()
