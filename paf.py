import re

re_tag = re.compile("(\S\S:[AZif]):(\S+)")
re_cs = re.compile("([:=*+-])(\d+|[A-Za-z]+)")

def get_edit_distance_pairs_from_tag(tag):
    seq = re.findall(re_cs, tag)
    pairs = []
    current_pos = 0
    for (operator, value) in seq:
        if operator == ':':
            for i in range(int(value)):
                pairs.append([current_pos, 0])
                current_pos = current_pos+1
        elif operator == '=':
            for i in range(len(value)):
                pairs.append([current_pos, 0])
                current_pos = current_pos+1
        elif operator == '-':
            for i in range(len(value)):
                pairs.append([current_pos, 1])
                current_pos = current_pos+1
        elif operator == '*':
            pairs.append([current_pos, 1])
            current_pos = current_pos+1
        elif operator == '+':
            for i in range(len(value)):
                if len(pairs) == 0:
                    pairs.append([current_pos, 0])
                pair = pairs.pop()
                pair[1] = pair[1]+1
                pairs.append(pair)
        else:
            raise ValueError("Unknown operator {}".format(operator))
    return pairs

def edit_distance_in_range(pairs, ref_start, ref_end):
    edit_distance = 0
    for pair in pairs:
        if pair[0] >= ref_start and pair[0] < ref_end:
            edit_distance = edit_distance + pair[1]
    return edit_distance

def get_pos_pairs_from_tag(tag, strand):
    seq = re.findall(re_cs, tag)
    pairs = []
    current_ref_pos = 0
    current_query_pos = 0
    step = 1
    if strand == "-":
        seq.reverse()
        step = -1
        current_ref_pos = calculate_ref_length(seq)-1
    for (operator, value) in seq:
        if operator == ':':
            for i in range(int(value)):
                pairs.append([current_ref_pos, current_query_pos])
                current_ref_pos = current_ref_pos+step
                current_query_pos = current_query_pos+1
        elif operator == '=':
            for i in range(len(value)):
                pairs.append([current_ref_pos, current_query_pos])
                current_ref_pos = current_ref_pos+step
                current_query_pos = current_query_pos+1
        elif operator == '-':
            for i in range(len(value)):
                pairs.append([current_ref_pos, current_query_pos])
                current_ref_pos = current_ref_pos+step
        elif operator == '*':
            pairs.append([current_ref_pos, current_query_pos])
            current_ref_pos = current_ref_pos+step
            current_query_pos = current_query_pos+1
        elif operator == '+':
            for i in range(len(value)):
                pairs.append([current_ref_pos, current_query_pos])
                current_query_pos = current_query_pos+1
        else:
            raise ValueError("Unknown operator: {}".format(operator))
    return pairs

def ref_pos_for_query_pos(pairs, query_pos):
    ref_pos_set = set()
    for pair in pairs:
        if pair[1] == query_pos:
            ref_pos_set.add(pair[0])
    return ref_pos_set

def query_pos_for_ref_pos(pairs, ref_pos):
    query_pos_set = set()
    for pair in pairs:
        if pair[0] == ref_pos:
            query_pos_set.add(pair[1])
    return query_pos_set

def calculate_ref_length(seq):
    length = 0
    for (operator, value) in seq:
        if operator == ':':
            length = length + int(value)
        elif operator == '=' or operator == '-':
            length = length + len(value)
        elif operator == '*':
            length = length + 1
        elif operator == '+':
            length = length + 0
        else:
            raise ValueError("Unknown operator {}".format(operator))
    return length

class Paf:
    def __init__(self, paf_line):
        fields = paf_line.split()
        self.qname = fields[0]
        self.qlength = int(fields[1])
        self.qstart = int(fields[2])
        self.qend = int(fields[3])
        self.strand = fields[4]
        self.rname = fields[5]
        self.rlength = int(fields[6])
        self.rstart = int(fields[7])
        self.rend = int(fields[8])
        self.numMatches = int(fields[9])
        self.alignmentLength = int(fields[10])
        self.mapQ = int(fields[11])
        self.tags = {}
        for i in range(12, len(fields)):
            m = re_tag.match(fields[i])
            self.tags[m.group(1)] = m.group(2)
        self.pairs = None
        self.edit_distance_pairs = None

    def get_tag(self, tag_name):
        try:
            tag = self.tags[tag_name]
        except:
            raise ValueError("{} tag not found".format(tag_name))
        return tag

    def get_pos_pairs(self):
        if self.pairs != None:
            return self.pairs
        tag = self.get_tag("cs:Z")
        self.pairs = get_pos_pairs_from_tag(tag, self.strand)
        return self.pairs

    def get_edit_distance_pairs(self):
        if self.edit_distance_pairs != None:
            return self.edit_distance_pairs
        tag = self.get_tag("cs:Z")
        self.edit_distance_pairs = get_edit_distance_pairs_from_tag(tag)
        return self.edit_distance_pairs

    def ref_pos_for_query_pos(self, query_pos):
        pairs = self.get_pos_pairs()
        return ref_pos_for_query_pos(pairs, query_pos)

    def query_pos_for_ref_pos(self, ref_pos):
        pairs = self.get_pos_pairs()
        return query_pos_for_ref_pos(pairs, ref_pos)

    def isInAlignedRangeRef(self, pos):
        return pos >= self.rstart and pos < self.rend

    def query_pos_for_ref_pos_in_paf(self, ref_pos):
        if ref_pos < self.rstart or ref_pos >= self.rend:
            raise ValueError("ref pos {} is not in aligned region: {}-{}".format(
                                ref_pos, self.rstart, self.rend))
        pos = ref_pos - self.rstart
        return_set = set()
        for query_pos in self.query_pos_for_ref_pos(pos):
            return_set.add(query_pos + self.qstart)
        return return_set

    def ref_pos_for_query_pos_in_paf_line(self, query_pos):
        if query_pos < self.qstart or query_pos >= self.qend:
            raise ValueError("query pos {} is not in aligned region: {}-{}".format(
                                query_pos, self.qstart, self.qend))
        pos = query_pos - self.qstart
        return_set = set()
        for ref_pos in self.ref_pos_for_query_pos(pos):
            return_set.add(ref_pos + self.rstart)
        return return_set

    def edit_distance_for_region_in_paf_line(self, ref_start, ref_end):
        if ref_start < self.rstart or ref_start >= self.rend:
            raise ValueError("start coordinate {} is not in aligned region: {}-{}".format(
                                ref_start, self.rstart, self.rend))
        if ref_end <= self.rstart or ref_end > self.rend:
            raise ValueError("end coordinate {} is not in aligned region: {}-{}".format(
                                ref_end, self.rstart, self.rend))
        if ref_start >= ref_end:
            raise ValueError("start coordinate is greater than end, please specify [start, end)")
        start = ref_start - self.rstart
        end = ref_end - self.rstart
        return edit_distance_in_range(self.get_edit_distance_pairs(), start, end)

