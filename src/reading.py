class ParsingDNA:
    all_lines = []

    def read_fasta_file(self, file):
        self.all_lines = [line.rstrip() for line in file]

    def fasta_to_dict_multi_line(self):
        """Преобразует fasta файл с многострочной последовательностью в словарь"""
        postions_name_seq = [i for i, elem in enumerate(self.all_lines) if elem.startswith('>')]
        seqs = {}
        for i, pos in enumerate(postions_name_seq):
            if i == len(postions_name_seq)-1:
                seqs[self.all_lines[pos][1:]] = ''.join(self.all_lines[pos+1:])
            else: seqs[self.all_lines[pos][1:]] = ''.join(self.all_lines[pos+1:postions_name_seq[i+1]])
        return seqs

    def fasta_to_dict_single_line(self):
        """Преобразует fasta файл с однострочной последовательностью в словарь"""
        name_seqs = [name_seq[1:] for i, name_seq in enumerate(self.all_lines) if i % 2 == 0]
        seqs = [seq for i, seq in enumerate(self.all_lines) if i % 2 != 0]
        seqs = {name_seqs[i]: seqs[i] for i in range(len(name_seqs))}
        return seqs

    def open_fasta_file(self, path_to_fasta_file: str, multi_line_format: bool = True):
        with open(path_to_fasta_file, 'r') as f:
            self.read_fasta_file(f)
            if multi_line_format:
                return self.fasta_to_dict_multi_line()
            else: return self.fasta_to_dict_single_line()

    def remove_coords_from_read_name_after_bedtools(self, reads: dict) -> dict:
        """
        Удаляет координаты транскрипта после bedtools getfasta
        >NCKAP5::chr2:133029574-133029594   =>     >NCKAP5
        """
        clean_keys = {}
        for key in reads:
            new_key = key.split("::")[0]
            clean_keys[new_key] = reads[key]
        return clean_keys