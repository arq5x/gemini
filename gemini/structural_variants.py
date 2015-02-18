class StructuralVariant(object):
    def __init__(self, var):
        self.var = var

    def is_precise(self):
        if self.var.INFO.get('IMPRECISE') is not None:
            return False
        else:
            return True

    def get_ci_left(self):
        cipos = self.var.INFO.get('CIPOS')
        if cipos is not None:
            ci_min, ci_max = [int(s) for s in cipos]
            return self.var.POS + ci_min, self.var.POS + ci_max
        else:
            return (None, None)

    def get_ci_right(self):
        ciend = self.var.INFO.get('CIEND')
        if ciend is not None:
            ci_min, ci_max = [int(s) for s in ciend]
            return self.var.end + ci_min, self.var.end + ci_max
        else:
            return (None, None)

    def get_sv_tool(self):
        return self.var.INFO.get('TOOL')

    def get_length(self):
        # for some reason, this is defined as a list
        length = self.var.INFO.get('SVLEN')
        if length is not None:
            # sometimes SVLEN is defined as a list in the VCF header,
            # yet othertime it is defined as an integer. handle both.
            if isinstance(length, list):
                return length[0]
            elif isinstance(length, int):
                return length

    def get_evidence_type(self):
        return self.var.INFO.get('EVTYPE')

    def get_event_id(self):
        return self.var.INFO.get('EVENT')

    def get_mate_id(self):
        return self.var.INFO.get('MATEID')

    def get_strand(self):
        """
        REF ALT Meaning
        s t[p[ piece extending to the right of p is joined after t
        s t]p] reverse comp piece extending left of p is joined after t
        s ]p]t piece extending to the left of p is joined before t
        s [p[t reverse comp piece extending right of p is joined before t

        """
        # multi-line SV
        if self.var.INFO.get('SVTYPE') == 'BND':
            if self.var.ALT[0][0] == '[': #[19:8195598[C
                return "-+"
            elif self.var.ALT[0][0] == ']': #]19:4529597]A
                return "--"
            elif self.var.ALT[0][-1] == '[': #A[19:8417020[
                return "++"
            elif self.var.ALT[0][-1] == ']': #T]19:8195491]
                return "+-"

        # single-line SV
        elif self.var.INFO.get('SVTYPE') == 'DEL':
            return "++"
        elif self.var.INFO.get('SVTYPE') == 'DUP':
            return "++"
        elif self.var.INFO.get('SVTYPE') == 'INV':
            return "+-"
        else:
            return None
