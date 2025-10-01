import os
from unittest import TestCase
from ddt import ddt, data, unpack

from src.bioomics import ParseImgtAnnot

@ddt
class TestParseImgtAnnot(TestCase):

    def setUp(self):
        self.data_dir = '/home/yuan/data/IMGT/3Dstructure-DB/IMGT3DFlatFiles/'
        self.obj = {}
        for inn in ('9447', '8127', '10636', '10562', '12254'):
            file_name = f'IMGT-MAB-{inn}.inn.gz'
            infile = os.path.join(self.data_dir, file_name)
            self.obj[inn] = ParseImgtAnnot(infile)

    @data(
        ['9447', True],
        ['8127', True],
        ['10636', True],
    )
    @unpack
    def test_env(self, inn_number, expect):
        file_name = f'IMGT-MAB-{inn_number}.inn.gz'
        infile = os.path.join(self.data_dir, file_name)
        result = os.path.isfile(infile)
        assert  result == expect

    @data(
        ['9447', 'INN number', ['9447']],
        ['8127', 'INN number', ['8127']],
        ['9447', 'Common name', ['IMC-RON8', 'RON8', 'IMC RON-8']],
        ['8127', 'Common name', None],
        ['9447', 'INN name', ['narnatumab']],
        ['8127', 'INN name', None],
        ['9447', 'IMGT receptor type', ['IG']],
        ['8127', 'IMGT receptor type', None],
        ['9447', 'IMGT receptor description', ['IG-GAMMA-1_KAPPA']],
        ['8127', 'IMGT receptor description', None],
        ['9447', 'Species', ['Homo sapiens (human)']],
        ['8127', 'Species', None],
        ['9447', 'Proposed list', ['L105 (2011)']],
        ['8127', 'Proposed list', ['L85 (2001)']],
        ['9447', 'Recommended list', ['R67 (2012)']],
        ['8127', 'Recommended list', ['R47 (2002)']],
        ['9447', 'Cas number', ['1188275-92-4_H 1188275-92-4_L']],
        ['8127', 'Cas number', None],
        ['9447', 'Molecular formula', ['C6454H10026N1754O2020S44']],
        ['8127', 'Molecular formula', None],
        ['9447', 'Glycosylation sites', ["302 302''"]],
        ['8127', 'Glycosylation sites', None],
        ['9447', 'Chain ID', ['9447_H', '9447_L']],
        ['8127', 'Chain ID', None],
        ['10636', 'Chain ID', ['10636_H', '10636_L', '10636_M', '10636_N']],
    )
    @unpack
    def test_next_line(self, inn, key, expect):
        c = self.obj[inn]
        result = c.next_line(key)
        assert result == expect

    @data(
        ['9447', 'Chain ID', '9447_H (9447H)'],
        ['9447', 'IMGT chain description', 'H-GAMMA-1'],
        ['9447', 'IMGT domain description', 'VH'],
        ['9447', 'wrong key', None],
        ['8127', 'Chain ID', None],
    )
    @unpack
    def test_line_2cols(self, inn, key, expect):
        c = self.obj[inn]
        result = c.line_2cols(key)
        assert result == expect


    @data(
        ['9447', '9447_H',
            'Chain ID                9447_H (9447H)',
            'C-DOMAIN      CSVMHEA.LHNHYTQKSLSLSP...'],
        ['9447', '9447_L',
            'Chain ID                9447_L (9447L)',
            'C-DOMAIN      CEVTHQG..LSSPVTKSFNRGEC',],
        ['9447', 'wrong id', None, None],
        ['10562', '10562_H',
            'Chain ID                10562_H (10562H)',
            'V-DOMAIN      QSYDSS..LSGWVFGGGTKLTVL',],
        ['8127', None, None, None],
    )
    @unpack
    def test_chain_block(self, inn, chain_id, expect_first, expect_last):
        c = self.obj[inn]
        result = c.chain_block(chain_id)
        if result:
            assert result[0] == expect_first
            assert result[-1] == expect_last
        else:
            assert result == []

    @data(
        ['9447', 'Intrachain disulfide bridges', {'H-H': \
            '22-96 149-205 266-326 372-430 22-96 149-205 266-326 372-430', \
                'L-L': '23-88 134-194 23-88 134-194'}],
        ['9447', 'Interchain disulfide bridges', {'H-H': \
            '234-234 231-231', 'H-L': '214-225 214-225', 'L-L': ''}],
        ['9447', 'wrong key', {}],
        ['8127', None, {}],
    )
    @unpack
    def test_block_lines(self, inn, key, expect):
        c = self.obj[inn]
        result = c.block_lines(key)
        assert result == expect

    def test_chain_vdomain_block(self):
        c = self.obj['9447']
        chain_block = c.chain_block('9447_H')
        result = c.chain_vdomain_block(chain_block)
        assert len(result) == 1
        assert result[0][0] == 'IMGT domain description  VH'
        assert result[0][10] == '                          [   CDR1   ]              '
        assert result[0][11] == 'EVQLVESGG.GLVQPGGSLRLSCAASGFTF....SSYLMTWVRQAPGKGLEW'
        assert len(result[0][10]) == len(result[0][11])
        assert result[0][0] == 'IMGT domain description  VH'
        assert result[0][12] == '   [  CDR2  ]                                       '
        assert result[0][13] == 'VANIKQD..GSEKYYVDSVK.GRFTISRDNAKNSLNLQMNSLRAEDTAVYYC'
        assert len(result[0][12]) == len(result[0][13])
        assert result[0][0] == 'IMGT domain description  VH'
        assert result[0][14] == '[    CDR3     ]           '
        assert result[0][15] == 'TRDGYSSGRHYGMDVWGQGTTVIVSS'
        assert len(result[0][14]) == len(result[0][15])

        chain_block = c.chain_block('9447_L')
        result = c.chain_vdomain_block(chain_block)
        assert len(result) == 1
        assert result[0][0] == 'IMGT domain description  V-KAPPA'
        assert result[0][-2] == '[   CDR3    ]          '
        assert result[0][-1] == 'QQRSN....WPRTFGQGTKVEIK'
        assert len(result[0][-2]) == len(result[0][-1])

        c = self.obj['10562']
        chain_block = c.chain_block('10562_H')
        result = c.chain_vdomain_block(chain_block)
        assert len(result) == 2
        assert result[0][0] == 'IMGT domain description  VH'
        assert result[0][-1] == 'AKMTSN.AFAFDYWGQGTLVTVSS'
        assert result[1][0] == 'IMGT domain description  V-LAMBDA'
        assert result[1][-1] == 'QSYDSS..LSGWVFGGGTKLTVL'

        c = self.obj['12254']
        chain_block = c.chain_block('12254_H')
        result = c.chain_vdomain_block(chain_block)
        assert len(result) == 2
        assert result[0][0] == 'IMGT domain description  V-KAPPA'
        assert result[0][-1] == 'QQRSS....YPYTFGGGTKVEIK'
        assert result[1][0] == 'IMGT domain description  VH'
        assert result[1][-1] == 'ATYGHYYGYMFAYWGQGTLVTVSS'

    def test_chain_cdomain_block(self):
        c = self.obj['9447']
        chain_block = c.chain_block('9447_H')
        result = c.chain_cdomain_block(chain_block)
        assert len(result) == 3
        assert result[0][0] == 'IMGT domain description  CH1'
        assert result[0][-1] == 'CNVNHKP..SNTKVDKRV'
        assert result[1][0] == 'IMGT domain description  CH2'
        assert result[1][-1] == 'CKVSNKA..LPAPIEKTISKAK'
        assert result[2][0] == 'IMGT domain description  CH3'
        assert result[2][-1] == 'CSVMHEA.LHNHYTQKSLSLSP...'
        chain_block = c.chain_block('9447_L')
        result = c.chain_cdomain_block(chain_block)
        assert len(result) == 1
        assert result[0][0] == 'IMGT domain description  C-KAPPA'
        assert result[0][-1] == 'CEVTHQG..LSSPVTKSFNRGEC'

        c = self.obj['10562']
        chain_block = c.chain_block('10562_H')
        result = c.chain_cdomain_block(chain_block)
        assert result == []


    @data(
        ['9447', 'EVQLVESGGGLVQPGGSLRLSCAASGFTF' + \
            'SSYLMTWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLNLQMN' + \
            'SLRAEDTAVYYCTRDGYSSGRHYGMDVWGQGTTVIVSSASTKGPSVFPLAPSSKS' + \
            'TSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVP' + \
            'SSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFP' + \
            'PKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNST' + \
            'YRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPS' + \
            'REEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSK' + \
            'LTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK'
        ],
        ['8127', ''],
    )
    @unpack
    def test_chain_seq(self, input, expect):
        c = self.obj[input]
        result = c.chain_sequence()
        assert result['chain_seq'] == expect
    
    @data(
        [
            '9447', [
            {'type': 'VH', 'start': 1, 'end': 122, 'seq': 'EVQLVESGGGLVQPGGSL' +\
                'RLSCAASGFTFSSYLMTWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSL' +\
                'NLQMNSLRAEDTAVYYCTRDGYSSGRHYGMDVWGQGTTVIVS', 'region': 'D1'},
            {'type': 'CH1', 'start': 123, 'end': 220, 'seq': 'ASTKGPSVFPLAPSS' +\
                'KSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSS' +\
                'LGTQTYICNVNHKPSNTKVDKR', 'region': 'D2'},
            {'type': 'HINGE-REGION'}, {'type': 'CH2', 'start': 236, 'end': 345, \
                'seq': 'APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVE' +\
                'VHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKA', 'region': 'D4'},
            {'type': 'CH3', 'start': 346, 'end': 450, 'seq': 'GQPREPQVYTLPPSR' +\
                'EEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKS' +\
                'RWQQGNVFSCSVMHEALHNHYTQKSLSLS', 'region': 'D5'}]
        ],
        ['8127', []],
    )
    @unpack
    def test_chain_regions(self, input, expect):
        c = self.obj[input]
        result = c.chain_sequence()
        assert result['regions'] == expect

    @data(
        ['9447', '9447_H', 1],
        ['9447', '9447_L', 1],
        ['10562', '10562_H', 2],
        ['12254', '12254_H', 2],
    )
    @unpack
    def test_vdomain_num(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        assert len(results) == expect

    @data(
        ['9447', '9447_H', [[
            {'name': 'CDR1', 'start': 26, 'end': 38, 'align': 'GFTF....SSYL', 'seq': 'GFTFSSYL'}, 
            {'name': 'CDR2', 'start': 55, 'end': 65, 'align': 'IKQD..GSEK', 'seq': 'IKQDGSEK'}, 
            {'name': 'CDR3', 'start': 104, 'end': 119, 'align': 'TRDGYSSGRHYGMDV', 'seq': 'TRDGYSSGRHYGMDV'}
        ]]],
        ['9447', '9447_L', [[
            {'name': 'CDR1', 'start': 26, 'end': 38, 'align': 'QSV......SRY', 'seq': 'QSVSRY'},
            {'name': 'CDR2', 'start': 55, 'end': 65, 'align': 'DA.......S', 'seq': 'DAS'},
            {'name': 'CDR3', 'start': 104, 'end': 117, 'align': 'QQRSN....WPRT', 'seq': 'QQRSNWPRT'}
        ]]],
        ['10562', '10562_H', [[
            {'name': 'CDR1', 'start': 26, 'end': 38, 'align': 'GFTF....RSYA', 'seq': 'GFTFRSYA'},
            {'name': 'CDR2', 'start': 55, 'end': 65, 'align': 'ISGR..GDNT', 'seq': 'ISGRGDNT'},
            {'name': 'CDR3', 'start': 104, 'end': 117, 'align': 'AKMTSN.AFAFDY', 'seq': 'AKMTSNAFAFDY'},
        ], [
            {'align': 'SSNIG...AGYG', 'end': 38, 'name': 'CDR1', 'seq': 'SSNIGAGYG', 'start': 26},
            {'align': 'GN.......T', 'end': 65, 'name': 'CDR2', 'seq': 'GNT', 'start': 55},
            {'align': 'QSYDSS..LSGWV', 'end': 117, 'name': 'CDR3', 'seq': 'QSYDSSLSGWV', 'start': 104}
        ]]],
        ['12254', '12254_H', [[
            {'name': 'CDR1', 'start': 26, 'end': 38, 'align': 'SSV.......SY', 'seq': 'SSVSY'},
            {'name': 'CDR2', 'start': 55, 'end': 65, 'align': 'ST.......S', 'seq': 'STS'},
            {'name': 'CDR3', 'start': 104, 'end': 117, 'align': 'QQRSS....YPYT', 'seq': 'QQRSSYPYT'}
        ], [
            {'name': 'CDR1', 'start': 26, 'end': 38, 'align': 'GYTF....TDYN', 'seq': 'GYTFTDYN'}, 
            {'name': 'CDR2', 'start': 55, 'end': 65, 'align': 'IYPN..NGGT', 'seq': 'IYPNNGGT'}, 
            {'name': 'CDR3', 'start': 104, 'end': 117, 'align': 'ATYGHYYGYMFAY', 'seq': 'ATYGHYYGYMFAY'}
        ]]],
    )
    @unpack
    def test_vdomain_cdr(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['cdr'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H', ['VH',]],
        ['9447', '9447_L', ['V-KAPPA',]],
        ['10562', '10562_H', ['VH', 'V-LAMBDA']],
        ['10562', '10562_L', []],
        ['12254', '12254_H', ['V-KAPPA', 'VH']],
    )
    @unpack
    def test_vdomain_desc(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['IMGT domain description'] for i in results]
        assert result == expect

    @data(
        [
            '9447',
            '9447_H', 
            [
                ['Homo sapiens IGHV3-7*01 (95.9%)', 'Homo sapiens IGHV3-7*02 (95.9%)',
                'Homo sapiens IGHV3-7*03 (95.9%)', 'Homo sapiens IGHJ6*01 (88.2%)',
                'Homo sapiens IGHJ6*02 (88.2%)',],
            ]
        ],
        [
            '9447',
            '9447_L',
            [
                ['Homo sapiens IGKV3-11*01 (98.9%)', 'Homo sapiens IGKJ1*01 (100%)'],
            ]
        ],
        [
            '10562',
            '10562_H',
            [
                ['Homo sapiens IGHV3-23*04 (94.9%)', 'Homo sapiens IGHJ4*01 (100%)',
                    'Homo sapiens IGHJ4*02 (100%)', 'Homo sapiens IGHJ4*03 (100%)'],
                ['Homo sapiens IGLV1-40*01 (97%)', 'Homo sapiens IGLJ3*02 (100%)']
            ]
        ],
        [
            '12254',
            '12254_H',
            [
                ['Mus musculus IGKV4-57*01 (87.1%)', 'Mus musculus IGKJ2*01 (91.7%)'],
                ['Mus musculus IGHV1S29*02 (81.4%)', 'Mus musculus IGHJ3*01 (92.9%)']
            ]
        ],
    )
    @unpack
    def test_vdomain_allele(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['IMGT gene and allele'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H', ['[8.8.15]',]],
        ['9447', '9447_L', ['[6.3.9]', ]],
        ['10562', '10562_H', ['[8.8.12]', '[9.3.11]']],
        ['12254', '12254_H', ['[5.3.9]', '[8.8.13]']],
    )
    @unpack
    def test_vdomain_cdr_length(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['CDR-IMGT length'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H', ['[25.17.38.11]',]],
        ['9447', '9447_L', ['[26.17.36.10]', ]],
        ['10562', '10562_H', ['[25.17.38.11]', '[25.17.36.10]']],
        ['12254', '12254_H', ['[25.17.36.10]', '[25.17.38.11]']],
    )
    @unpack
    def test_vdomain_fr_length(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['FR-IMGT length'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYLMTWVRQAPGKGLEWVAN' + \
                    'IKQDGSEKYYVDSVKGRFTISRDNAKNSLNLQMNSLRAED' + \
                    'TAVYYCTRDGYSSGRHYGMDVWGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L', 
            [
                'EIVLTQSPATLSLSPGERATLSCRASQSVSRYLAWYQQKPGQAPRLL' + \
                    'IYDASNRATGIPARFSGSGSGTDFTLTISSLEPE' + \
                    'DFAVYYCQQRSNWPRTFGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGGGLVQPGGSLRLSCAASGFTFRSYAMSWVRQAPGKGLEWVSAISGR' + \
                    'GDNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKMTSN' + \
                    'AFAFDYWGQGTLVTVSS',
                'QSVLTQPPSVSGAPGQRVTISCTGSSSNIGAGYGVHWYQQLPGTAPKLL' + \
                    'IYGNTNRPSGVPDRFSGFKSGTSASLAITGLQAEDEADYYCQSYDSS' + \
                    'LSGWVFGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                'IQLTQSPSSLSASPGDRVTITCSASSSVSYMHWFQQKPGKAPKLWIYSTSNLA' + \
                    'SGVPARFSGSGSGTSYSLTISRLQPEDIATYYCQQRSSYPYTFGGGTKVEIK', 
                'EVQLVQSGAEVKKPGASVKVSCKASGYTFTDYNMDWVKQSPGQGLEWMGYIYPN' + \
                    'NGGTGYNQKFKSKVTITVDTSTSTAYMELHSLRSEDTAVYYCATYGHY' + \
                    'YGYMFAYWGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGG.GLVQPGGSLRLSCAASGFTF' + \
                    '....SSYLMTWVRQAPGKGLEWVANIKQD..GSEKYYVDSVK.GRFTISRDNAKNSLNLQ' + \
                    'MNSLRAEDTAVYYCTRDGYSSGRHYGMDVWGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L',
            [
                'EIVLTQSPATLSLSPGERATLSCRASQSV' + \
                    '......SRYLAWYQQKPGQAPRLLIYDA.......SNRATGIP.ARFSGSG..SGTDF' + \
                    'TLTISSLEPEDFAVYYCQQRSN....WPRTFGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGG.GLVQPGGSLRLSCAASGFTF....' + \
                    'RSYAMSWVRQAPGKGLEWVSAISGR..GDNTYYADSVK.GRFTISRDNSKNTLYLQMNSLRAEDTA' + \
                    'VYYCAKMTSN.AFAFDYWGQGTLVTVSS',
                'QSVLTQPPS.VSGAPGQRVTISCTGSSSNIG' + \
                    '...AGYGVHWYQQLPGTAPKLLIYGN.......TNRPSGVP.DRFSGFK..SGTSASLAITGL' + \
                    'QAEDEADYYCQSYDSS..LSGWVFGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                '.IQLTQSPSSLSASPGDRVTITCSASSSV.......' + \
                    'SYMHWFQQKPGKAPKLWIYST.......SNLASGVP.ARFSGSG..SGTSYSLTISRLQPEDIATYYC' + \
                    'QQRSS....YPYTFGGGTKVEIK', 
                'EVQLVQSGA.EVKKPGASVKVSCKASGYTF....' + \
                    'TDYNMDWVKQSPGQGLEWMGYIYPN..NGGTGYNQKFK.SKVTITVDTSTSTAYMELHSLRSEDTA' + \
                    'VYYCATYGHYYGYMFAYWGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq_align(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq_align'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGGGLVQPGGSLRLSCAAS[GFTFSSYL]MTWVRQAPGKGLEWVAN[IKQDGSEK]YYVDSVKGRFTISRD' + \
                    'NAKNSLNLQMNSLRAEDTAVYYC[TRDGYSSGRHYGMDV]WGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L',
            [
                'EIVLTQSPATLSLSPGERATLSCRAS[QSVSRY]LAWYQQKPGQAPRLLIY[DAS]NRATGIPARFSGSG' + \
                    'SGTDFTLTISSLEPEDFAVYYC[QQRSNWPRT]FGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGGGLVQPGGSLRLSCAAS[GFTFRSYA]MSWVRQAPGKGLEWVSA[ISGRGDNT]YYADSVKGRFTISRDNSK' + \
                    'NTLYLQMNSLRAEDTAVYYC[AKMTSNAFAFDY]WGQGTLVTVSS',
                'QSVLTQPPSVSGAPGQRVTISCTGS' + \
                    '[SSNIGAGYG]VHWYQQLPGTAPKLLIY[GNT]NRPSGVPDRFSGFK' + \
                    'SGTSASLAITGLQAEDEADYYC[QSYDSSLSGWV]FGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                'IQLTQSPSSLSASPGDRVTITCSAS[SSVSY]MHWFQQKPGKAPKLWIY[STS]NLASGVPARFSGSGSG' + \
                    'TSYSLTISRLQPEDIATYYC[QQRSSYPYT]FGGGTKVEIK', 
                'EVQLVQSGAEVKKPGASVKVSCKAS[GYTFTDYN]MDWVKQSPGQGLEWMGY[IYPNNGGT]GYNQKFKSKVTITVDTST' + \
                    'STAYMELHSLRSEDTAVYYC[ATYGHYYGYMFAY]WGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq_cdr(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq_cdr'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGG.GLVQPGGSLRLSCAAS' + \
                    '[GFTF....SSYL]MTWVRQAPGKGLEWVAN[IKQD..GSEK]YYVDSVK.GRFTISRD' + \
                    'NAKNSLNLQMNSLRAEDTAVYYC[TRDGYSSGRHYGMDV]WGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L',
            [
                'EIVLTQSPATLSLSPGERATLSCRAS' + \
                    '[QSV......SRY]LAWYQQKPGQAPRLLIY[DA.......S]NRATGIP.ARFSGSG.' + \
                    '.SGTDFTLTISSLEPEDFAVYYC[QQRSN....WPRT]FGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGG.GLVQPGGSLRLSCAAS' + \
                    '[GFTF....RSYA]MSWVRQAPGKGLEWVSA[ISGR..GDNT]YYADSVK.GRFTISRDNSK' + \
                    'NTLYLQMNSLRAEDTAVYYC[AKMTSN.AFAFDY]WGQGTLVTVSS',
                'QSVLTQPPS.VSGAPGQRVTISCTGS' + \
                    '[SSNIG...AGYG]VHWYQQLPGTAPKLLIY[GN.......T]NRPSGVP.DRFSGFK.' + \
                    '.SGTSASLAITGLQAEDEADYYC[QSYDSS..LSGWV]FGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                '.IQLTQSPSSLSASPGDRVTITCSAS' + \
                    '[SSV.......SY]MHWFQQKPGKAPKLWIY[ST.......S]NLASGVP.ARFSGSG..SG' + \
                    'TSYSLTISRLQPEDIATYYC[QQRSS....YPYT]FGGGTKVEIK', 
                'EVQLVQSGA.EVKKPGASVKVSCKAS' + \
                    '[GYTF....TDYN]MDWVKQSPGQGLEWMGY[IYPN..NGGT]GYNQKFK.SKVTITVDTST' + \
                    'STAYMELHSLRSEDTAVYYC[ATYGHYYGYMFAY]WGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq_align_cdr(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq_align_cdr'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGG.GLVQPGGSLRLSCAAS<CDR1>MTWVRQAPGKGLEWVAN<CDR2>YYVDSVK.GRFTISRD' + \
                    'NAKNSLNLQMNSLRAEDTAVYYC<CDR3>WGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L',
            [
                'EIVLTQSPATLSLSPGERATLSCRAS<CDR1>LAWYQQKPGQAPRLLIY<CDR2>NRATGIP.ARFSGSG.' + \
                    '.SGTDFTLTISSLEPEDFAVYYC<CDR3>FGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGG.GLVQPGGSLRLSCAAS<CDR1>MSWVRQAPGKGLEWVSA<CDR2>YYADSVK.GRFTISRDNSK' + \
                    'NTLYLQMNSLRAEDTAVYYC<CDR3>WGQGTLVTVSS',
                'QSVLTQPPS.VSGAPGQRVTISCTGS<CDR1>VHWYQQLPGTAPKLLIY<CDR2>NRPSGVP.DRFSGFK.' + \
                    '.SGTSASLAITGLQAEDEADYYC<CDR3>FGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                '.IQLTQSPSSLSASPGDRVTITCSAS<CDR1>MHWFQQKPGKAPKLWIY<CDR2>NLASGVP.ARFSGSG..SG' + \
                    'TSYSLTISRLQPEDIATYYC<CDR3>FGGGTKVEIK', 
                'EVQLVQSGA.EVKKPGASVKVSCKAS<CDR1>MDWVKQSPGQGLEWMGY<CDR2>GYNQKFK.SKVTITVDTST' + \
                    'STAYMELHSLRSEDTAVYYC<CDR3>WGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq_align_cdr_mask(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq_align_cdr_mask'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'EVQLVESGGGLVQPGGSLRLSCAAS<CDR1>MTWVRQAPGKGLEWVAN<CDR2>YYVDSVKGRFTISRD' + \
                    'NAKNSLNLQMNSLRAEDTAVYYC<CDR3>WGQGTTVIVSS',
            ]
        ],
        ['9447', '9447_L',
            [
                'EIVLTQSPATLSLSPGERATLSCRAS<CDR1>LAWYQQKPGQAPRLLIY<CDR2>NRATGIPARFSGSG' + \
                    'SGTDFTLTISSLEPEDFAVYYC<CDR3>FGQGTKVEIK',
            ]
        ],
        ['10562', '10562_H',
            [
                'QVQLVESGGGLVQPGGSLRLSCAAS<CDR1>MSWVRQAPGKGLEWVSA<CDR2>YYADSVKGRFTISRDNSK' + \
                    'NTLYLQMNSLRAEDTAVYYC<CDR3>WGQGTLVTVSS',
                'QSVLTQPPSVSGAPGQRVTISCTGS<CDR1>VHWYQQLPGTAPKLLIY<CDR2>NRPSGVPDRFSGFK' + \
                    'SGTSASLAITGLQAEDEADYYC<CDR3>FGGGTKLTVL'
            ]
        ],
        ['12254', '12254_H',
            [
                'IQLTQSPSSLSASPGDRVTITCSAS<CDR1>MHWFQQKPGKAPKLWIY<CDR2>NLASGVPARFSGSGSG' + \
                    'TSYSLTISRLQPEDIATYYC<CDR3>FGGGTKVEIK', 
                'EVQLVQSGAEVKKPGASVKVSCKAS<CDR1>MDWVKQSPGQGLEWMGY<CDR2>GYNQKFKSKVTITVDTST' + \
                    'STAYMELHSLRSEDTAVYYC<CDR3>WGQGTLVTVSS'
            ]
        ],
    )
    @unpack
    def test_vdomain_seq_cdr_mask(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        results = c.chain_v_domain(vdomain_blocks)
        result = [i['seq_cdr_mask'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H', ['CH1', 'CH2', 'CH3',]],
        ['9447', '9447_L', ['C-KAPPA',]],
    )
    @unpack
    def test_cdomain_desc(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_cdomain_block(chain_block)
        results = c.chain_c_domain(vdomain_blocks)
        result = [i['IMGT domain description'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H', 
            [
                ['Homo sapiens IGHG1*01 (99%)', 'Homo sapiens IGHG1*02 (99%)',
                    'Homo sapiens IGHG1*03 (100%)',],
                ['Homo sapiens IGHG1*01 (100%)', 'Homo sapiens IGHG1*02 (100%)',
                    'Homo sapiens IGHG1*03 (100%)',],
                ['Homo sapiens IGHG1*01 (98.1%)', 'Homo sapiens IGHG1*02 (98.1%)',
                    'Homo sapiens IGHG1*03 (100%)',]
            ]
        ],
        ['9447', '9447_L', [['Homo sapiens IGKC*01 (100%)',]]]
    )
    @unpack
    def test_cdomain_allele(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_cdomain_block(chain_block)
        results = c.chain_c_domain(vdomain_blocks)
        result = [i['IMGT gene and allele'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                '....ASTKGPSVFPLAPSSKSTS...' + \
                    'GGTAALGCLVKDYFP..EPVTVSWNSGALTS....GVHTFPAVLQSS......GLYSLSSVVTV' + \
                    'PSSSL...GTQTYICNVNHKP..SNTKVDKRV',
                '..APELLGGPSVFLFPPKPKDTLMI.SRTPE' + \
                    'VTCVVVDVSHEDPEVKFNWYVDGVEVH...NAKTKPREEQYN......STYRVVSVLTV' + \
                    'LHQDW..LNGKEYKCKVSNKA..LPAPIEKTISKAK',
                '....GQPREPQVYTLPPSREEMT...KNQ' + \
                    'VSLTCLVKGFYP..SDIAVEWESNGQPEN...NYKTTPPVLDSD......GSFFLYSKLTV' + \
                    'DKSRW..QQGNVFSCSVMHEA.LHNHYTQKSLSLSP...',
            ]
        ],
        ['9447', '9447_L',
            [
                '....RTVAAPSVFIFPPSDEQLK...SGTASV' + \
                    'VCLLNNFYP..REAKVQWKVDNALQSG..NSQESVTEQDSKD.....STYSLSSTLTL' + \
                    'SKADY..EKHKVYACEVTHQG..LSSPVTKSFNRGEC',
            ]
        ],
    )
    @unpack
    def test_cdomain_seq_align(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_cdomain_block(chain_block)
        results = c.chain_c_domain(vdomain_blocks)
        result = [i['seq_align'] for i in results]
        assert result == expect

    @data(
        ['9447', '9447_H',
            [
                'ASTKGPSVFPLAPSSKSTS' + \
                    'GGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTV' + \
                    'PSSSLGTQTYICNVNHKPSNTKVDKRV',
                'APELLGGPSVFLFPPKPKDTLMISRTPE' + \
                    'VTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTV' + \
                    'LHQDWLNGKEYKCKVSNKALPAPIEKTISKAK',
                'GQPREPQVYTLPPSREEMTKNQVSLTCLVK' + \
                    'GFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHE' + \
                    'ALHNHYTQKSLSLSP',

            ]
        ],
        ['9447', '9447_L',
            [
                'RTVAAPSVFIFPPSDEQLKSGTASV' + \
                    'VCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTL' + \
                    'SKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC'
            ]
        ],
    )
    @unpack
    def test_cdomain_seq(self, input_id, input_chain, expect):
        c = self.obj[input_id]
        chain_block = c.chain_block(input_chain)
        vdomain_blocks = c.chain_cdomain_block(chain_block)
        results = c.chain_c_domain(vdomain_blocks)
        result = [i['seq'] for i in results]
        assert result == expect



    def test_v_domain_empty(self):
        c = self.obj['8127']
        chain_block = c.chain_block('wrong')
        assert chain_block == []
        vdomain_blocks = c.chain_vdomain_block(chain_block)
        assert vdomain_blocks == []
        result = c.chain_v_domain(vdomain_blocks)
        assert result == []

    def test_c_domain_empty(self):
        c = self.obj['8127']
        result = c.chain_c_domain([])
        assert result == []

        c = self.obj['10562']
        chain_block = c.chain_block('10562_H')
        cdomain_blocks = c.chain_cdomain_block(chain_block)
        assert cdomain_blocks == []
        result = c.chain_c_domain(cdomain_blocks)
        assert result == []