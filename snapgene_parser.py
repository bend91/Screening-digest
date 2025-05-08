import struct


class SnapGeneObject:
    def __init__(self, filepath):
        self._filepath = filepath
        self._meta_data = {}
        self._dna = {}
        self._notes = ""
        self.parse_file(filepath)
        ...

    def parse_file(self):
        fileobject = open(self._filepath, "rb")
        if fileobject.read(1) != b"\t":
            print("Not a snapgene file")

        length = unpack(fileobject, 4, "I")
        title = fileobject.read(8).decode("ascii")
        if title != "SnapGene":
            print("Not a snapgene file")
        self._meta_data["isDNA"] = unpack(fileobject, 2, "H")
        self._meta_data["exportVersion"] = unpack(fileobject, 2, "H")
        self._meta_data["importVersion"] = unpack(fileobject, 2, "H")
        self._meta_data["other"] = []
        parsing = True
        while parsing:
            next_byte = fileobject.read(1)
            if next_byte == b"":
                print("End of file")
                parsing = False
                continue
            block_size = unpack(fileobject, 4, "I")
            if ord(next_byte) == 0:
                props = unpack(fileobject, 1, "b")
                self._dna["topology"] = "circular" if props & 0x01 else "linear"
                self._dna["strandedness"] = "double" if props & 0x02 > 0 else "single"
                self._dna["damMethylated"] = props & 0x04 > 0
                self._dna["dcmMethylated"] = props & 0x08 > 0
                self._dna["ecoKIMethylated"] = props & 0x010 > 0
                self._dna["length"] = block_size - 1
                self._dna["sequence"] = fileobject.read(block_size - 1).decode("ascii")
            elif ord(next_byte) == 6:
                block_content = fileobject.read(block_size).decode("utf-8")
                note_data = parse_dict(block_content)
                self._notes = note_parser(note_data)
            else:
                self._meta_data["other"].append(fileobject.read(block_size))


def note_parser(notes):
    parsed_note = {}
    note_offset = 0
    note_sub = notes[notes.find(">") + 1:notes.find("</Notes")]
    while len(note_sub) > 1:
        header = note_sub[note_sub.find("<")+1:note_sub.find(">")]
        note_sub = note_sub[note_sub.find(">")+1:]
        content = note_sub[:note_sub.find("<")]
        note_sub = note_sub[note_sub.find(">")+1:]
        parsed_note[header] = content
    return parsed_note


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
    extra_data = []
    parsing = True
    while parsing:
        next_byte = fileobject.read(1)
        if next_byte == b"":
            print("End of file")
            parsing = False
            continue
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
            data["notes"] = note_parser(note_data)
        else:
            # props = unpack(fileobject, 1, "b")
            # extra_data.append(props)
            extra_data.append(fileobject.read(block_size))

fileobject.close()
    return data
