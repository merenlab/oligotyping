# v.110811

class SequenceSource:
    def __init__(self, f_name):
        self.pos = 0
        self.id  = None
        self.seq = None
        self.file_pointer = open(f_name)
        self.total_seq = len([l for l in self.file_pointer.readlines() if l.startswith('>')])
        self.file_pointer.seek(0)

    def next(self):
        self.id = self.file_pointer.readline()[1:].strip()
        self.seq = None

        sequence = ''
        while 1:
            line = self.file_pointer.readline()
            if not line:
                if len(sequence):
                    self.seq = sequence
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            sequence += line.strip()

        self.seq = sequence
        self.pos += 1
        return True

    def reset(self):
        self.pos = 0
        self.id  = None
        self.seq = None
        self.file_pointer.seek(0)



