DEFAULT_PIPE =\
{
'wget' : {
'key':'wget',
'description':'Download command',
'switch':True,
'path':r'',
'outdir':'fastq',
'output_ext':'',
'inputdir':'',
'input_ext':'',
'structure':''
},
'trim' : {
'key':'trim_galore',
'description':'Trims adapters',
'switch':True,
'path':r'',
'outdir':'post_trim',
'output_ext':'.fq.gz',
'inputdir':'fastq',
'input_ext':'.fastq',
'structure':'key settings fastq_pre'
},
'fastqc':{
'key':'fastqc',
'description':'Create a readable quality report in html format',
'switch':True,
'path':r'fastqc',
'outdir':'post_trim',
'output_ext':'',
'inputdir':'',
'input_ext':'',
'structure':''
},
'fastqc_pre':{
'key':'fastqc',
'switch':True,
'description':'Create a readable quality report in html format',
'path':r'fastqc',
'outdir':'fastq',
'output_ext':'',
'inputdir':'',
'input_ext':'',
'structure':''
},
'bwa-align' : {
'key':'bwa mem',
'description':'Align reads to reference',
'switch':True,
'path':r'',
'outdir':'alignment',
'output_ext':'.sam',
'inputdir':'post_trim',
'input_ext':'.fq.gz',
'structure':'key settings fasta fastq output_pipe'
},
'view':{
'key':'samtools view',
'description':'Convert sam to bam',
'switch':True,
'path':r'',
'outdir':'alignment',
'output_ext':'.bam',
'inputdir':'alignment',
'input_ext':'.sam',
'structure':'key settings _sam output_pipe'
},
'index':{
'key':'samtools index',
'description':'Index bam for fast access via samtools',
'switch':True,
'path':r'',
'outdir':'alignment',
'output_ext':'.bam.bai',
'inputdir':'alignment',
'input_ext':'.bam',
'structure':''
},
'mpileup':{
'key':'samtools mpileup',
'description':'Create a readable mpileup format',
'switch':True,
'path':r'',
'outdir':'pileup',
'output_ext':'.pileup',
'inputdir':'alignment',
'input_ext':'.bam',
'structure':'key settings command_input output_pipe'
},
'sort':{
'key': 'samtools sort',
'description':'Sort reads',
'switch':True,
'path':r'',
'outdir':'',
'output_ext':'',
'inputdir':'alignment',
'input_ext':'.bam',
'structure':'key settings command_input'
},
'htseq-count':{
'key': 'htseq-count',
'description': 'Measure expression',
'switch': False,
'path': r'',
'outdir':'count',
'output_ext':'.txt',
'inputdir':'alignment',
'input_ext':'.bam',
'structure': 'key settings input gtf output_pipe'
},
'STAR-alignment':{
'key':'STAR',
'description':'Align reads to ref',
'switch': False,
'path': 'STAR',
'outdir':'',
'output_ext':'',
'inputdir':'post_trim',
'input_ext':'.fq.gz',
'structure':'key settings'
}
}

DEFAULT_SETTINGS =\
{
'trim':{
'quality':{'description':'Trim low-quality ends from reads', 'key':'-q', 'value':20, 'switch':True},
'outdir':{'description':'Choose output directory', 'key':'-o', 'value':'_fastq_dir', 'switch':True}
},
'bwa-align':{
'seed':{'description':'Minimum seed length', 'key':'-k', 'value':19, 'switch':True},
},
'view':{
'header':{'description':'Include the header in the output', 'key':'-h', 'value':'', 'switch':True},
'bam':{'description':'Output in BAM format', 'key':'-b', 'value':'', 'switch':True}
},
'index':{},
'mpileup':{
'quality':{'description':'Minimum base quality for a base to be considered', 'key':'-Q', 'value':30, 'switch':True},
'reference':{'description':'Fasta type ref file','key':'-f', 'value':'_fasta', 'switch':True}
},
'fastqc':{
'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':'_fastq_dir', 'switch':True}
},
'fastqc_pre':{
'outdir':{'description':'Output directory', 'key' : '--outdir', 'value':f'_fastq_pre_dir', 'switch':True}
},
'wget':{},
'sort':{
'output':{'description':'Output directory', 'key' : '-o', 'value': '_bam', 'switch':True}
},
'htseq-count':{
'ID':{'description': 'ID attribute to use', 'key' : '-i', 'value' : 'ID', 'switch':True},
'Additional attribute':{'description': 'Additional attribute to use', 'key':'--additional-attr', 'value':'gene','switch':True},
'Feature column name':{'description': 'Feature column name to use','key':'-t', 'value':'gene', 'switch':True},
'stranded':{'description': 'Whether the reads are stranded', 'key':'-s','value':'no','switch':True},
'format':{'description':'Format to read input in','key':'-f','value':'bam','switch':True},
'ambiguous':{'description':'How to handle ambiguous reads', 'key':'--nonunique','value':'all','switch':True}
},
'STAR-alignment':{
'genomeDir':
{
'description':'Input ref',
'key' : '--genomeDir',
'value':'_fasta_dir',
'switch':True
},
'readFilesIn':
{
'description':'Input fastq',
'key': '--readFilesIn',
'value':'_fastq',
'switch':True
},
"outSAMtype":
{
"description": "output format",
"key": "--outSAMtype",
"switch": True,
"value": "BAM SortedByCoordinate"
},
"readFilesCommand":
{
"description": "Input format",
"key": "--readFilesCommand",
"switch": True,
"value": "gunzip -c"
},
"runThreadN":
{
"description": "The amount of threads to use",
"key": "--runThreadN",
"switch": True,
"value": 4
},
"outFileNamePrefix":
{
"description": "Output prefix",
"key": '--outFileNamePrefix',
'switch': True,
'value': '_STAR_out'
},
'limit sort RAM':
{
'description':'whatever',
'key':'--limitBAMsortRAM',
'value':1205130643,
'switch':True
}}
}
DEFAULT_ORDER = ['wget','fastqc_pre','trim','fastqc','bwa-align','view','sort','index','mpileup']

DEFAULT_MAIN = \
  {
  'ignore_paired_end': [False, 'Ignore paired-end reads'],
  'qsub_queue' : ['bioinfo.q', 'choose which cluster queue to send the jobs into'],
  'environment':['bash', 'choose the linux environment to use'],
  'free_form':['', 'any command you type here will be added to every single .csh file created, enter nothing to wipe setting']
  }