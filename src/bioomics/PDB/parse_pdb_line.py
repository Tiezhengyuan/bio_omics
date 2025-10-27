
import re

class ParsePdbLine:
    def __init__(self, line:str):
        self.line = line
    
    def parse_atom(self):
        rec = {
            'record_type': self.line[:3],
            'serial_number': self.line[6:11],
            'atom_name': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
            '': self.line[:],
        }


13-16	Atom name	left*	character
17	Alternate location indicator		character
18-20§	Residue name	right	character
22	Chain identifier		character
23-26	Residue sequence number	right	integer
27	Code for insertions of residues		character
31-38	X orthogonal Å coordinate	right	real (8.3)
39-46	Y orthogonal Å coordinate	right	real (8.3)
47-54	Z orthogonal Å coordinate	right	real (8.3)
55-60	Occupancy	right	real (6.2)
61-66	Temperature factor	right	real (6.2)
73-76	Segment identifier¶	left	character
77-78	Element symbol	right	character
79-80	Charge