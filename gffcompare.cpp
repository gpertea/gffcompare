#include "GArgs.h"
#include <ctype.h>
#include <errno.h>
#include "gtf_tracking.h"

#define VERSION "0.12.6"

#define USAGE "Usage:\n\
gffcompare [-r <reference_mrna.gtf> [-R]] [-T] [-V] [-s <seq_path>]\n\
    [-o <outprefix>] [-p <cprefix>] \n\
    {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}\n\
\n\
 GffCompare provides classification and reference annotation mapping and\n\
 matching statistics for RNA-Seq assemblies (transfrags) or other generic\n\
 GFF/GTF files.\n\
 GffCompare also clusters and tracks transcripts across multiple GFF/GTF\n\
 files (samples), writing matching transcripts (identical intron chains) into\n\
 <outprefix>.tracking, and a GTF file <outprefix>.combined.gtf which \n\
 contains a nonredundant set of transcripts across all input files (with\n\
 a single representative transfrag chosen for each clique of matching transfrags\n\
 across samples).\n\
\n\
 Options:\n\
 -v display gffcompare version (also --version)\n\
 -i provide a text file with a list of (query) GTF files to process instead\n\
    of expecting them as command line arguments (useful when a large number\n\
    of GTF files should be processed)\n\
\n\
 -r reference annotation file (GTF/GFF)\n\
 -R for -r option, consider only the reference transcripts that\n\
    overlap any of the input transfrags (Sn correction)\n\
 -Q for -r option, consider only the input transcripts that\n\
    overlap any of the reference transcripts (Precision correction);\n\
    (Warning: this will discard all \"novel\" loci!)\n\
 -M discard (ignore) single-exon transfrags and reference transcripts\n\
 -N discard (ignore) single-exon reference transcripts\n\
 -D discard \"duplicate\" query transfrags (i.e. same intron chain) within\n\
    a single sample (disable \"annotation\" mode for a single file); \n\
    this option is automatically enabled when multiple query files are provided\n\
 -S when -D is enabled (or multiple query files are provided), perform a more \n\
    strict duplicate checking: only discard matching (same intron chain) query \n\
    transcripts from the same sample if their boundaries are fully contained \n\
    within (or same with) matching transcripts\n\
    if --strict-match is also given, exact match of all exons is required\n\
 --no-merge : disable close-exon merging (default: merge exons separated by\n\
   \"introns\" shorter than 5 bases\n\
\n\
 -s path to genome sequences (optional); this can be either a multi-FASTA\n\
    file or a directory containing single-fasta files (one for each contig);\n\
    repeats must be soft-masked (lower case) in order to be able to classify\n\
    transfrags as repeats\n\
\n\
 -e when estimating exon level accuracy, this is the maximum range\n\
    variation allowed for the free ends of terminal exons (default 100);\n\
	this terminal exon restriction  applies to transcript level accuracy\n\
	when --strict-match option is given\n\
 --strict-match : transcript matching takes into account the -e range\n\
    for terminal exons; code '=' is only assigned if transcript ends are\n\
    within that range, otherwise code '~' is assigned just for intron chain\n\
    match (or significant overlap in the case of single exon transcripts)\n\
\n\
 -d max. distance (range) for grouping transcript start sites (100)\n\
 -V verbose processing mode (also shows GFF parser warnings)\n\
 -T do not generate .tmap and .refmap files for each input file\n\
 --chr-stats: the .stats file will show summary and accuracy data\n\
   per reference contig/chromosome\n\
 -j <output.tab> if -r was given, writes novel junctions in this file\n\
 --debug : enables -V and generates additional files: \n\
    <outprefix>.Q_discarded.lst, <outprefix>.missed_introns.gff,\n\
    <outprefix>.R_missed.lst\n\
\n\
Options for the combined GTF output file:\n\
 -p the name prefix to use for consensus transcripts in the \n\
    <outprefix>.combined.gtf file (default: 'TCONS')\n\
 -C discard matching and \"contained\" transfrags in the GTF output\n\
    (i.e. collapse intron-redundant transfrags across all query files)\n\
 -A like -C but does not discard intron-redundant transfrags if they start\n\
    with a different 5' exon (keep alternate TSS)\n\
 -X like -C but also discard contained transfrags if transfrag ends stick out\n\
    within the container's introns\n\
 -K for -C/-A/-X, do NOT discard any redundant transfrag matching a reference\n\
"
/*
 --gids     : append related reference gene_id values to gene_id\n\
 --gidnames : append related reference gene_name values to gene_id\n\
     Note: --gids and --gidnames options are mutually exclusive!\n\
 --gnames   : make gene_name include reference gene_name values, if any\n\
*/

bool perContigStats=false; // -S to enable separate stats for every contig/chromosome
//bool generic_GFF=false;
//true if -G: won't discard intron-redundant transfrags
bool discardContained=false; //activated by -C, -A or -X
//bool showContained=true; // opposite of -C, default is to show them in combined.gtf
//-- moved to gtf_tracking
bool keepAltTSS=false; //-A option
bool keepRefMatching=false; //-K with -C/-A/-X
bool allowIntronSticking=false; //-X option
bool reduceQrys=false; //-Q
bool checkFasta=false;
bool tmapFiles=true;
//ref gene_id or gene_name values are separated by '|' (pipe character) when
//appended to the original gene_id
//gene_name values are separated by ',' (comma) in the gene_name attribute

//TODO implement these
bool set_gene_name=false; //gene_name set to the list of overlapping ref gene_names
bool gid_add_ref_gids=false; //append overlapping ref gene_ids to gene_id
bool gid_add_ref_gnames=false; //append overlapping ref gene_names to gene_id

bool qDupDiscard=false;

//strictMatching=false; // really match *all* exon coords for '=' class code!
          // use '~' class code for intron-chain matches
          // (or fuzzy single-exon transcripts overlaps)
//noMergeCloseExons=false; //prevent joining close exons?
bool only_spliced_refs=false;
int debugCounter=0;
bool gffAnnotate=false;
GStr consGTF;
int outConsCount=0;
int polyrun_range=2000; //polymerase run range 2KB
double scoreThreshold=0;
char* cprefix=NULL;
FILE* ffasta=NULL; //genomic seq file
FILE *f_ref=NULL; //reference mRNA GFF, if provided
FILE* f_in=NULL; //sequentially, each input GFF file
FILE* f_out=NULL; //stdout if not provided
GFastaHandler gfasta;
int xlocnum=0;
int tsscl_num=0; //for tss cluster IDs
int protcl_num=0; //for "unique" protein IDs within TSS clusters
int tssDist=100;
uint exonEndRange=100; // -e value, only used for exon level Sn/Sp
//int total_tcons=0;
int total_xloci_alt=0;

void openfwrite(FILE* &f, GArgs& args, char opt) {
  GStr s=args.getOpt(opt);
  if (!s.is_empty()) {
      if (s=='-')
       f=stdout;
      else {
       f=fopen(s,"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

//-- structure to keep track of data from multiple qry input files for a single genomic seq
class GSeqTrack {
  int gseq_id;
 public:
  const char* gseq_name;
  GList<GLocus>* rloci_f; //reference loci for this genomic sequence
  GList<GLocus>* rloci_r;
  GList<GXLocus> xloci_f; // extended super-loci across all qry datasets
  GList<GXLocus> xloci_r; // extended super-loci across all qry datasets
  GList<GXLocus> xloci_u; // extended super-loci across all qry datasets
  GVec<GSeqData*> qdata; //fixed order array with GSeqData for each qry input
                 //element in array is NULL if a qry file has no transcripts on this genomic sequence
  int get_gseqid() { return gseq_id; }

  GSeqTrack(int numqryfiles, int gid):xloci_f(true,true,false),
        xloci_r(true,true,false), xloci_u(true,true,false), qdata() {
    gseq_id=gid;
    if (gseq_id>=0) {
      gseq_name=GffObj::names->gseqs.getName(gseq_id);
      }
    rloci_f=NULL;
    rloci_r=NULL;
    //for (int i=0;i<MAX_QFILES;i++) qdata[i]=NULL;
    qdata.Resize(numqryfiles, NULL);
    }

  bool operator==(GSeqTrack& d){
      return (gseq_id==d.gseq_id);
      }
  bool operator<(GSeqTrack& d){
     return (gseq_id<d.gseq_id);
     }
};

//char* getFastaFile(int gseq_id);

// ref globals
bool haveRefs=false;  //true if a reference annotation (-r) is provide

GList<GSeqData> ref_data(true,true,true); //list of reference mRNAs and loci data for each genomic seq
              //each locus will keep track of any superloci which includes it, formed during the analysis

void processLoci(GSeqData& seqdata, GSeqData* refdata=NULL, int qfidx=0);

void reportStats(FILE* fout, const char* setname, GSuperLocus& stotal,
       GSeqData* seqdata=NULL, GSeqData* refdata=NULL, int qfidx=-1);

GSeqData* getQryData(int gid, GList<GSeqData>& qdata);
void trackGData(int qcount, GList<GSeqTrack>& gtracks, GStr& fbasename, FILE** ftr, FILE** frs);

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)
#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)

FILE* f_mintr=NULL; //missed ref introns (debug only)
FILE* f_qdisc=NULL; //discarded query transfrags (debug only)
FILE* f_rmiss=NULL; //missed ref transcripts (debug only)
FILE* f_nj=NULL; //-j novel junctions output file

bool multiexon_only=false;
bool multiexonrefs_only=false;

GHash<GStr*> refdescr;
void loadRefDescr(const char* fname);

GList<GStr> qryfiles(false,true,false);

//list of GSeqTrack data, sorted by gseq_id
GList<GSeqTrack> gseqtracks(true,true,true);
GSeqTrack* findGSeqTrack(int gsid);


int cmpGTrackByName(const pointer p1, const pointer p2) {
 return strcmp(((GSeqTrack*)p1)->gseq_name, ((GSeqTrack*)p2)->gseq_name);
}


void show_version() {
  GMessage("gffcompare v%s\n", VERSION);
}

void show_usage() {
  GMessage("gffcompare v%s\n", VERSION);
  GMessage( "-----------------------------\n");
  GMessage("%s\n", USAGE);
  }

void RefReqCheck(bool v, const char* opt) {
	if (v)
		GError("%s\nError: option %s requires reference annotation (-r)\n",
			USAGE, opt);
}

//preallocate the cmpovl slots for ref loci,
// to keep the info about the query loci overlapping each of them by qfidx
void prepRefLoci(GList<GLocus> & rloci, int qcount) {
	for (int i=0;i<rloci.Count();i++)
      for (int q=0;q<qcount;q++)
    	  rloci[i]->qlocovls.Add(new GList<GLocus>(true, false, true));
}


int main(int argc, char* argv[]) {
#ifdef HEAPROFILE
  if (!IsHeapProfilerRunning())
      HeapProfilerStart("./gffcompare_dbg.hprof");
#endif

  GArgs args(argc, argv,
		  "version;help;debug;gids;gidnames;gnames;no-merge;strict-match;"
		  "chr-stats;vACDSGEFJKLMNQTVRXhp:e:d:s:i:j:n:r:o:");
  int e;
  if ((e=args.isError())>0) {
    show_usage();
    GMessage("Invalid argument: %s\n", argv[e]);
    exit(1);
    }
  if (args.getOpt('h') || args.getOpt("help")){
    show_usage();
    exit(0);
    }
  if (args.getOpt('v') || args.getOpt("version")){
    show_version();
    exit(0);
  }
  debug=(args.getOpt("debug")!=NULL);
  tmapFiles=(args.getOpt('T')==NULL);
  multiexon_only=(args.getOpt('M')!=NULL);
  multiexonrefs_only=(args.getOpt('N')!=NULL);
  set_gene_name=(args.getOpt("gnames")!=NULL);
  gid_add_ref_gids=(args.getOpt("gids")!=NULL);
  gid_add_ref_gnames=(args.getOpt("gidnames")!=NULL);
  qDupDiscard=(args.getOpt('D')!=NULL);
  qDupStrict=(args.getOpt('S')!=NULL);
  stricterMatching=(args.getOpt("strict-match")!=NULL);
  //if (stricterMatching) exonEndRange=0;
  noMergeCloseExons=(args.getOpt("no-merge")!=NULL);
  if (gid_add_ref_gids && gid_add_ref_gnames)
	GError("Error: options --gids and --gidnames are mutually exclusive!\n");
  perContigStats=(args.getOpt("chr-stats")!=NULL);
  checkFasta=(args.getOpt('J')!=NULL);
  gtf_tracking_verbose=((args.getOpt('V')!=NULL) || debug);
  FILE* finlst=NULL;
  GStr s=args.getOpt('i');
  if (!s.is_empty()) {
	  if (s=='-')
		  finlst=stdin;
	  else {
		  finlst=fopen(s,"r");
		  if (finlst==NULL) GError("Error opening file: %s\n", s.chars());
	  }
  }

  if (finlst) {
	  GLineReader* lr=new GLineReader(finlst);
	  char* l=NULL;
	  while ((l=lr->getLine())!=NULL) {
		  if (strlen(l)<2 || startsWith(l,"# ") || isspace(*l)) continue;
		  if (!fileExists(l)) GError("Error: cannot locate input file: %s\n", l);
		  qryfiles.Add(new GStr(l));
	  }
	  delete lr;
	  //if (qryfiles.Count()>10)
	  gtf_tracking_largeScale=true;
	  tmapFiles=false;
  }
  else {
	  numQryFiles=args.startNonOpt();
	  char *infile=NULL;
	  if (numQryFiles>0) {
		  if (numQryFiles>6) {
			  gtf_tracking_largeScale=true;
			  tmapFiles=false;
		  }
		  while ((infile=args.nextNonOpt())!=NULL) {
			  if (!fileExists(infile)) GError("Error: cannot locate input file: %s\n", infile);
			  qryfiles.Add(new GStr(infile));
		  } //for each argument
	  }
  }
  numQryFiles=qryfiles.Count();
  if (numQryFiles==0) {
	  show_usage();
	  exit(1);
  }
  gfasta.init(args.getOpt('s'));
   // determine if -s points to a multi-fasta file or a directory
  //s=args.getOpt('c');
  //if (!s.is_empty()) scoreThreshold=s.asReal();
  s=args.getOpt('p');
  if (!s.is_empty()) cprefix=Gstrdup(s.chars());
               else  cprefix=Gstrdup("TCONS");
  s=args.getOpt('e');
  if (!s.is_empty())
	  exonEndRange=s.asInt();
  if (stricterMatching)
	  terminalMatchRange=exonEndRange;

  s=args.getOpt('d');
  if (!s.is_empty()) tssDist=s.asInt();

  s=args.getOpt('n');
  if (!s.is_empty()) loadRefDescr(s.chars());
  reduceRefs=(args.getOpt('R')!=NULL);
  reduceQrys=(args.getOpt('Q')!=NULL);

  //if a full pathname is given
  //the other common output files will still be created in the current directory:
  // .loci, .tracking, .stats
  GStr outbasename; //include path, if provided
  GStr outprefix; //without path and/or extension
  GStr outstats=args.getOpt('o');
  if (outstats.is_empty() || outstats=="-") {
       outstats="gffcmp";
       }
  outbasename=outstats;
  GStr outext(getFileExt(outstats.chars()));
  if (outext.is_empty()) {
    outext="stats";
    outstats.append(".stats");
    outbasename=outstats;
    }
    else outext.lower();
  if (outext=="txt" || outext=="out" || outext=="stats" || outext=="summary") {
      outbasename.cut(outbasename.length()-outext.length()-1);
  }

  outprefix=outbasename;
  int di=outprefix.rindex(CHPATHSEP);
  if (di>=0) outprefix.cut(0,di+1);
  if (gtf_tracking_verbose) GMessage("Prefix for output files: %s\n", outprefix.chars());

  s=args.getOpt('r');
  if (s.is_empty()) {
	RefReqCheck(multiexonrefs_only, "-M");
	RefReqCheck(reduceRefs, "-R");
	RefReqCheck(reduceQrys, "-Q");
	RefReqCheck(set_gene_name, "--gnames");
	RefReqCheck(gid_add_ref_gids, "--gids");
	RefReqCheck(gid_add_ref_gnames, "--gidnames");
  } else { //reference annotation file provided
    f_ref=fopen(s,"r");
    if (f_ref==NULL) GError("Error opening reference gff: %s\n",s.chars());
    haveRefs=true;
    if (gtf_tracking_verbose) GMessage("Loading reference transcripts..\n");
    read_mRNAs(f_ref, ref_data, &ref_data, true, -1, s.chars(),
    		(multiexonrefs_only || multiexon_only));
    haveRefs=(ref_data.Count()>0);
    //if (gtf_tracking_verbose) GMessage("..reference annotation loaded\n");
    if (haveRefs && !gtf_tracking_largeScale) {
    	//prepare qlocovls for each ref locus
    	for (int i=0;i<ref_data.Count();i++) {
    		prepRefLoci(ref_data[i]->loci_f, qryfiles.Count());
    		prepRefLoci(ref_data[i]->loci_r, qryfiles.Count());
    	}
    }

  }
  if (args.getOpt('C')) discardContained=true;
  if (args.getOpt('A')) {
	  //redundancy check will NOT consider containment for alt. TSS
	  keepAltTSS=true;
	  if (discardContained)
		   GMessage("Warning: option -C will be ignored, -A takes precedence.\n");
	  discardContained=true;
  }
  if (args.getOpt('X')) {
	discardContained=true;
	allowIntronSticking=true;
  }
  keepRefMatching=(args.getOpt('K')!=NULL);
  if (keepRefMatching && !discardContained) {
	  GMessage("Warning: -K option ignored, requires -C, -A or -X\n");
	  keepRefMatching=false;
  }
  s=args.getOpt('j');
  if (!s.is_empty()) {
     f_nj=fopen(s.chars(),"w");
     if (f_nj==NULL) GError("Error creating file %s!\n",s.chars());
  }
  if (debug) { //create a few more files potentially useful for debugging
        s=outbasename;
        s.append(".missed_introns.gff");
        f_mintr=fopen(s.chars(),"w");
        if (f_mintr==NULL) GError("Error creating file %s!\n",s.chars());

        s=outbasename;
        s.append(".R_missed.gff");
        f_rmiss=fopen(s.chars(),"w");
        if (f_rmiss==NULL) GError("Error creating file %s!\n",s.chars());

        if (reduceQrys) { //only if -Q option was used
           s=outbasename;
           s.append(".Q_discarded.lst");
           f_qdisc=fopen(s.chars(),"w");
           if (f_qdisc==NULL) GError("Error creating file %s!\n",s.chars());
        }
  }

  f_out=fopen(outstats, "w");
  if (f_out==NULL) GError("Error creating output file %s!\n", outstats.chars());
  fprintf(f_out, "# gffcompare v%s | Command line was:\n#", VERSION);
  for (int i=0;i<argc-1;i++)
    fprintf(f_out, "%s ", argv[i]);
  fprintf(f_out, "%s\n#\n", argv[argc-1]);
  //int qfileno=0;
  GList<GSeqData>** qrysdata=NULL;
  FILE** tfiles=NULL;
  FILE** rtfiles=NULL;
  GMALLOC(qrysdata, numQryFiles*sizeof(GList<GSeqData>*));
  if (tmapFiles) {
	  GMALLOC(tfiles, numQryFiles*sizeof(FILE*));
	  if (haveRefs) {
		  GMALLOC(rtfiles, numQryFiles*sizeof(FILE*));
	  }
  }
  gffAnnotate=(numQryFiles==1 && !discardContained && haveRefs && !qDupDiscard && !qDupStrict);
  consGTF=outbasename;
  if (gffAnnotate) consGTF.append(".annotated.gtf");
              else consGTF.append(".combined.gtf");

  for (int fi=0;fi<qryfiles.Count();fi++) {
    GStr in_file(qryfiles[fi]->chars());
    GStr infname(getFileName(qryfiles[fi]->chars())); //file name only
    GStr indir(qryfiles[fi]->chars());
    di=indir.rindex(CHPATHSEP);
    if (di>=0) indir.cut(di+1); //directory path for this input file
          else indir=""; //current directory

    //if (debug || (gtf_tracking_verbose && !gtf_tracking_largeScale))
    if (qryfiles.Count()>1)
        GMessage("Loading query file #%d: %s\n",fi+1, in_file.chars());
    if (in_file=="-") { f_in=stdin; in_file="stdin"; }
      else {
        f_in=fopen(in_file.chars(),"r");
        if (f_in==NULL)
            GError("Cannot open input file %s!\n",in_file.chars());
        }
    //f_in is the query gff file to process
    GStr sbase(indir);
    sbase.append(outprefix);
    sbase.append(".");
    sbase.append(infname);
    if (tmapFiles) {
        //-- we should keep the infname path, otherwise the remaining file names
        //   may be the same and clobber each other
        s=sbase;
        s.append(".tmap");
        tfiles[fi]=fopen(s.chars(),"w");
        if (tfiles[fi]==NULL)
          GError("Error creating file '%s'!\n",s.chars());
        fprintf(tfiles[fi],"ref_gene_id\tref_id\tclass_code\tqry_gene_id\tqry_id\tnum_exons\tFPKM\tTPM\tcov\tlen\tmajor_iso_id\tref_match_len\n");
        if (haveRefs) {
          s=sbase;
          s.append(".refmap");
          rtfiles[fi]=fopen(s.chars(),"w");
          if (rtfiles[fi]==NULL)
             GError("Error creating file '%s'!\n",s.chars());
          fprintf(rtfiles[fi],"ref_gene\tref_id\tclass_code\tqry_id_list\n");
        }
    }
    GList<GSeqData>* pdata=new GList<GSeqData>(true,true,true);
    qrysdata[fi]=pdata;
    //int discard_check=discard_redundant;
    //if (keepRefMatching) {
    //  discard_check=0;
    //}
    read_mRNAs(f_in, *pdata, &ref_data, !gffAnnotate, fi, in_file.chars(), multiexon_only);
    GSuperLocus gstats;
    for (int g=0;g<pdata->Count();g++) { //for each genomic sequence in this qry dataset
        int gsid=pdata->Get(g)->get_gseqid();
        GSeqData* refdata=getRefData(gsid, ref_data);//ref data for this contig
        if (!gtf_tracking_largeScale)
          processLoci(*(pdata->Get(g)), refdata, fi);
        GSeqTrack* seqtrack=findGSeqTrack(gsid); //this will add a gseqtrack if it doesn't exist
        // for gsid
        if (refdata!=NULL) {
          seqtrack->rloci_f=&(refdata->loci_f);
          seqtrack->rloci_r=&(refdata->loci_r);
        }
        seqtrack->qdata[fi]=pdata->Get(g);
        //will only gather data into stats if perContig==false
        if (!gtf_tracking_largeScale) reportStats(f_out, getGSeqName(gsid), gstats,
              pdata->Get(g), refdata, fi);
    } //for each genomic sequence data for the current query file
    //-- there could also be genomic sequences with no qry transcripts
    //   but only reference transcripts so they weren't found above
    if (haveRefs && !reduceRefs && !gtf_tracking_largeScale) {
        for (int r=0;r<ref_data.Count();r++) {
          GSeqData* refdata=ref_data[r];
          int gsid=refdata->get_gseqid();
          if (getQryData(gsid, *pdata)==NULL) {
            reportStats(f_out, getGSeqName(gsid), gstats, NULL, refdata);
            }//completely missed all refdata on this contig
        }
    }
    //now report the summary:
    if (!gtf_tracking_largeScale) reportStats(f_out, in_file.chars(), gstats);
      //qfileno++;
  }//for each query file
  if (f_mintr!=NULL) fclose(f_mintr); //to write missed introns
  if (f_qdisc!=NULL) fclose(f_qdisc); //to write discarded query transcripts
  if (f_rmiss!=NULL) fclose(f_rmiss); //to write missed references
  if (f_nj!=NULL) fclose(f_nj); //to write missed introns
  gseqtracks.setSorted(&cmpGTrackByName);
  if (gtf_tracking_verbose && numQryFiles>1)
	   GMessage("Tracking transcripts across %d query file(s)..\n", numQryFiles);
  trackGData(numQryFiles, gseqtracks, outbasename, tfiles, rtfiles);
  fprintf(f_out, "\n Total union super-loci across all input datasets: %d \n", xlocnum);
  if (numQryFiles>1) {
      fprintf(f_out, "  (%d multi-transcript, ~%.1f transcripts per locus)\n",
           total_xloci_alt, ((double)(GXConsensus::count))/xlocnum);
      }
  int redundant_consensi=GXConsensus::count-outConsCount;
  fprintf(f_out, "%d out of %d consensus transcripts written in %s (%d discarded as redundant)\n",
		  outConsCount, GXConsensus::count, consGTF.chars(), redundant_consensi);
  if (gtf_tracking_verbose) GMessage("Cleaning up..\n");
  GFREE(cprefix);
  // clean up
  for (int i=0;i<numQryFiles;i++) {
    delete qrysdata[i];
  }
  GFREE(qrysdata);
  GFREE(tfiles);
  GFREE(rtfiles);
  gseqtracks.Clear();
  FWCLOSE(f_out);
  if (gtf_tracking_verbose) GMessage("Done.\n");
  ref_data.Clear();
  //getchar();
} //main ends here

void show_exons(FILE* f, GffObj& m) {
  fprintf(f,"(");
  int imax=m.exons.Count()-1;
  for (int i=0;i<=imax;i++) {
    if (i==imax) fprintf(f,"%d-%d)",m.exons[i]->start, m.exons[i]->end);
            else fprintf(f,"%d-%d,",m.exons[i]->start, m.exons[i]->end);
    }
}

bool exon_match(GXSeg& r, GXSeg& q, uint fuzz=0) {
 uint sd = (r.start>q.start) ? r.start-q.start : q.start-r.start;
 uint ed = (r.end>q.end) ? r.end-q.end : q.end-r.end;
 uint ex_range=exonEndRange;
 if (ex_range<=fuzz) ex_range=fuzz;
 if ((r.flags&1) && (q.flags&1)) { // first exon ?
	if (sd>ex_range) return false;
 }
 else {
	if (sd>fuzz) return false;
 }
 if ((r.flags&2) && (q.flags&2)) { // last exon ?
	if (ed>ex_range) return false;
 }
 else {
	if (ed>fuzz) return false;
 }
 return true;
}

void addQLocOvl(GLocus* rloc, GLocus* qloc, int qfidx) {
	//adds a qloc to a list of qloc overlaps (cmpovl list) based on qfidx
	//only applies to a ref locus
	if (qfidx>=rloc->qlocovls.Count())
		GError("Error: addQLocOvl() not ready!\n");
	rloc->qlocovls[qfidx]->Add(qloc);
}

void compareLoci2R(GList<GLocus>& loci, GList<GSuperLocus>& cmpdata,
                             GList<GLocus>& refloci, int qfidx) {
 cmpdata.Clear();//a new list of superloci will be built
 if (refloci.Count()==0 || loci.Count()==0) return;
 //reset cmpovl and stats
 for (int i=0;i<refloci.Count();i++)
	 refloci[i]->creset();
 //find loci with overlapping refloci
 //and store cmpovl links both ways for ALL loci and refloci on this strand
 for (int l=0;l<loci.Count();l++) {
   GLocus* locus=loci[l];
   locus->creset();
   for (int j=0;j<refloci.Count();j++) {
     if (refloci[j]->start>locus->end) {
         if (refloci[j]->start > locus->end + GFF_MAX_LOCUS) break;
         continue;
     }
     if (locus->start>refloci[j]->end) continue;
     //must check for proper exon overlap:
     if (locus->exonOverlap(*refloci[j])) {
    	 //cmpovl is a list of exon-overlapping loci
        locus->cmpovl.Add(refloci[j]); //qry locus adds this overlapping ref locus
        refloci[j]->cmpovl.Add(locus); //ref locus adds this overlapping qry locus
        addQLocOvl(refloci[j], locus, qfidx);
     }
   }//for each reflocus
 } //for each locus

 //create corresponding "superloci" from transitive overlapping between loci and ref
 for (int l=0;l<loci.Count();l++) {
  if (loci[l]->v!=0) continue; //skip, already processed
  GSuperLocus* super=new GSuperLocus();
  super->qfidx=qfidx;
  //try to find all other loci connected to this locus loci[l]
  GPVec<GLocus> lstack(false);  //traversal stack
  lstack.Push(loci[l]);
  while (lstack.Count()>0) {
      GLocus* locus=lstack.Pop();
      if (locus->v!=0 || locus->cmpovl.Count()==0) continue;
      super->addQlocus(*locus);
      locus->v=1;
      for (int r=0;r<locus->cmpovl.Count();r++) {
        GLocus* rloc=locus->cmpovl[r];
        if (rloc->v==0) {
          super->addRlocus(*rloc);
          rloc->v=1;
          for (int ll=0;ll<rloc->cmpovl.Count();ll++) {
              if (rloc->cmpovl[ll]->v==0)
            	  lstack.Push(rloc->cmpovl[ll]);
          }
        }
      } //for each overlapping ref locus
  } //while linking

  if (super->qloci.Count()==0) {
    delete super;
    continue; //try next query locus
  }
  //--here we have a "superlocus" region data on both qry and ref
  // -- analyze mexons matching (base level metrics)
  cmpdata.Add(super);
  //make each ref locus keep track of all superloci containing it
  for (int x=0;x<super->rmexons.Count();x++) {
    super->rbases_all += super->rmexons[x].end-super->rmexons[x].start+1;
  }
  for (int x=0;x<super->qmexons.Count();x++) {
    super->qbases_all += super->qmexons[x].end-super->qmexons[x].start+1;
  }
  int i=0; //locus mexons
  int j=0; //refmexons
  while (i<super->qmexons.Count() && j<super->rmexons.Count()) {
     uint istart=super->qmexons[i].start;
     uint iend=super->qmexons[i].end;
     uint jstart=super->rmexons[j].start;
     uint jend=super->rmexons[j].end;
     if (iend<jstart) { i++; continue; }
     if (jend<istart) { j++; continue; }
     //v--overlap here:
     uint ovlstart = jstart>istart? jstart : istart;
     uint ovlend = iend<jend ? iend : jend;
     uint ovlen=ovlend-ovlstart+1;
     super->baseTP+=ovlen; //qbases_cov
     if (iend<jend) i++;
               else j++;
  } //while mexons ovl search
  /* if (reduceRefs) {
    super->baseFP=super->qbases_all-super->baseTP;
    super->baseFN=super->rbases_all-super->baseTP;
    }
  */
  // -- exon level comparison:
  int* qexovl; //flags for qry exons with ref overlap
  GCALLOC(qexovl,super->quexons.Count()*sizeof(int));
  int* rexovl; //flags for ref exons with qry overlap
  GCALLOC(rexovl,super->ruexons.Count()*sizeof(int));
  for (int i=0;i<super->quexons.Count();i++) {
	  uint istart=super->quexons[i].start;
	  uint iend=super->quexons[i].end;
	  for (int j=0;j<super->ruexons.Count();j++) {
		  uint jstart=super->ruexons[j].start;
		  uint jend=super->ruexons[j].end;
		  if (iend<jstart) break;
		  if (jend<istart) continue;
		  //--- overlap here between quexons[i] and ruexons[j]
		  qexovl[i]++;
		  rexovl[j]++;
		  /*if (exon_match(super->quexons[i], super->ruexons[j],5)) {
			  if (!ATPfound) { //count a ref approx match only once
				  super->exonATP++;
				  ATPfound=true;
			  }*/
			  if (exon_match(super->quexons[i], super->ruexons[j])) {
				  if ((super->ruexons[j].flags & 4)==0) {
				      super->exonTP++;
				      super->ruexons[j].flags |= 4;
				  }
				  if ((super->quexons[i].flags & 4)==0) {
				      super->exonQTP++;
				      super->quexons[i].flags |= 4;
				  }
			  } //exact match
		  //} //fuzzy match
	  } //ref uexon loop
  } //qry uexon loop
  //DEBUG only:
  //if (super->exonTP > super->total_rexons)
  //   GMessage("Warning: superlocus %d-%d has exonTP %d > total rexons %d\n",
  //      super->start, super->end, super->exonTP, super->total_rexons);
  super->m_exons=0; //ref exons with no query overlap
  super->w_exons=0; //qry exons with no ref overlap
  for (int x=0;x<super->quexons.Count();x++)
       if (qexovl[x]==0) super->w_exons++;
  for (int x=0;x<super->ruexons.Count();x++)
       if (rexovl[x]==0) super->m_exons++;
  GFREE(rexovl);
  GFREE(qexovl);

  //-- intron level stats:
  //query:
  int* qinovl=NULL; //flags for qry introns with at some ref intron overlap
  int* qtpinovl=NULL; //flags for qry introns with ref intron match
  if (super->qintrons.Count()>0) {
    GCALLOC(qinovl,super->qintrons.Count()*sizeof(int));
    GCALLOC(qtpinovl,super->qintrons.Count()*sizeof(int));
  }
  //-- reference:
  int* rinovl=NULL; //flags for ref introns with qry overlap
  int* rtpinovl=NULL; //ref introns with qry intron match
  if (super->rintrons.Count()>0) {
    GCALLOC(rinovl,super->rintrons.Count()*sizeof(int));
    GCALLOC(rtpinovl,super->rintrons.Count()*sizeof(int));
  }
  for (int i=0;i<super->qintrons.Count();i++) {
    uint istart=super->qintrons[i].start;
    uint iend=super->qintrons[i].end;
    for (int j=0;j<super->rintrons.Count();j++) {
      uint jstart=super->rintrons[j].start;
      uint jend=super->rintrons[j].end;
      if (iend<jstart) break;
      if (jend<istart) continue;
      //--- overlap here between qintrons[i] and rintrons[j]
      qinovl[i]++;
      rinovl[j]++;
      if (super->qintrons[i].coordMatch(&super->rintrons[j])) {
         super->intronTP++;
         qtpinovl[i]++;
         rtpinovl[j]++;
      } //exact match
    } //ref intron loop
  } //qry intron loop
  super->m_introns=0; //ref introns with no query overlap (missed introns)
  super->w_introns=0; //qry introns with no ref overlap (wrong introns)
  for (int x=0;x<super->qintrons.Count();x++) {
       if (qinovl[x]==0) {
    	   super->w_introns++;
           //qry introns with no ref intron overlap at all
           super->i_qwrong.Add(super->qintrons[x]);
       } else if (qtpinovl[x]==0) { //novel introns = qry introns with no ref intron match
             super->i_qnotp.Add(super->qintrons[x]);
       }
  }
  for (int x=0;x<super->rintrons.Count();x++) {
       if (rinovl[x]==0) { //no qry intron overlap at all
             super->m_introns++; //ref introns totally missed = not having any query intron overlaps
             super->i_missed.Add(super->rintrons[x]);
       } else if (rtpinovl[x]==0) { //no qry intron match
    	    //ref introns with with no qry intron match
            super->i_notp.Add(super->rintrons[x]);
       }
  }
  GFREE(rinovl);
  GFREE(rtpinovl);
  GFREE(qinovl);
  GFREE(qtpinovl);

  // ---- now intron-chain and transcript matching
  GVec<char> matched_refs(super->rmrnas.Count(), '\0'); //keep track of matched refs
                              //'=' when "exact" (within terminalMatchRange), or if no strict matching was requested
                              //'~' when stricter transcript matching is activated and only the intron chain was matched
  //GVec<int> amatched_refs(super->rmrnas.Count(), 0); //keep track of fuzzy-matched refs
  for (int i=0;i<super->qmrnas.Count();i++) {
	  uint istart=super->qmrnas[i]->exons.First()->start;
	  uint iend=super->qmrnas[i]->exons.Last()->end;
	  for (int j=0;j<super->rmrnas.Count();j++) {
		  if (matched_refs[j]=='=') continue; //already counted as ichainTP and mrnaTP
		  uint jstart=super->rmrnas[j]->exons.First()->start;
		  uint jend=super->rmrnas[j]->exons.Last()->end;
		  if (iend<jstart) break;
		  if (jend<istart) continue;
		  //--- overlapping  transcripts ---
		  if (super->qmrnas[i]->udata & 2) continue; //already found a matching ref for this
		  GLocus* qlocus=((CTData*)super->qmrnas[i]->uptr)->locus;
		  GLocus* rlocus=((CTData*)super->rmrnas[j]->uptr)->locus;
		  int ovlen=0;
		  //look for a transcript match ('=' code for full exact exons match, '~' )
		  char tmatch=transcriptMatch(*(super->qmrnas[i]),*(super->rmrnas[j]), ovlen, terminalMatchRange);
		  //bool isTMatch=(tmatch>0);
		  if (tmatch) {
			  if (!stricterMatching) tmatch='=';
			  //at least the intron chains match !
			  if (super->qmrnas[i]->exons.Count()>1) {
				  super->ichainTP++;
				  qlocus->ichainTP++;
				  if ((super->qmrnas[i]->udata & 4) == 0) {
					  super->qmrnas[i]->udata |= 4;
				  }
				  if (matched_refs[j]==0) {
					  rlocus->ichainTP++;
					  matched_refs[j]=tmatch;
				  }
			  }
			  if (tmatch=='=') { //"full" or strict match
				  super->mrnaTP++;
				  qlocus->mrnaTP++;
				  rlocus->mrnaTP++;
				  if ((super->qmrnas[i]->udata & 2) ==0) {
					  super->qmrnas[i]->udata|=2;
				  }
				  matched_refs[j]='=';
			  }
		  } // ichain match found
	  } //ref loop
  } //qry loop
  //-- show missed references (not "matched" either '~' or '=') if requested
  if (f_rmiss!=NULL) {
	  for (int i=0;i<super->rmrnas.Count();i++) {
		  if (matched_refs[i]==0)
			  super->rmrnas[i]->printGxf(f_rmiss, pgffAny);
	  }
  }
  for (int ql=0;ql<super->qloci.Count();ql++) {
      if (super->qloci[ql]->mrnaTP>0)
                 super->locusQTP++;
  }
  for (int rl=0;rl<super->rloci.Count();rl++) {
      if (super->rloci[rl]->mrnaTP >0 )
                 super->locusTP++;
  }
 }//for each unlinked locus
}

//look for qry data for a specific genomic sequence
GSeqData* getQryData(int gid, GList<GSeqData>& qdata) {
  int qi=-1;
  GSeqData f(gid);
  GSeqData* q=NULL;
  if (qdata.Found(&f,qi))
        q=qdata[qi];
  return q;
}

const char* findDescr(GffObj* gfobj) {
  if (refdescr.Count()==0) return NULL;
  GStr* s=refdescr.Find(gfobj->getID());
  if (s==NULL) {
       s=refdescr.Find(gfobj->getGeneName());
       if (s==NULL) s=refdescr.Find(gfobj->getGeneID());
       }
  if (s!=NULL)
     return s->chars();
  return NULL;
}

const char* getGeneID(GffObj* gfobj) {
 //returns anything that might resemble a gene identifier for this transcript
 //or, if everything fails, returns the transcript ID
 const char* s=gfobj->getGeneID();
 if (s) return s;
 if ((s=gfobj->getGeneName())!=NULL) return s;
 if ((s=gfobj->getAttr("Name"))!=NULL) return s;
 return gfobj->getID();
}

const char* getGeneNameID(GffObj& gfobj) {
 //returns anything that might resemble a gene name or gene identifier for the transcript
 //or, if everything fails, returns the transcript ID
 const char* s=gfobj.getGeneName();
 if (s) return s;
 if ((s=gfobj.getGeneID())!=NULL) return s;
 if ((s=gfobj.getAttr("Name"))!=NULL) return s;
 return gfobj.getID();
}

const char* getGeneID(GffObj& gfobj) {
 return getGeneID(&gfobj);
}

void writeLoci(FILE* f, GList<GLocus> & loci) {
 for (int l=0;l<loci.Count();l++) {
   GLocus& loc=*(loci[l]);
   fprintf(f,"%s\t%s[%c]%d-%d\t", loc.mrna_maxcov->getID(),
       loc.mrna_maxcov->getGSeqName(),
           loc.mrna_maxcov->strand, loc.start,loc.end);
   //now print all transcripts in this locus, comma delimited
   int printfd=0;
   for (int i=0;i<loc.mrnas.Count();i++) {
      if (loc.mrnas[i]==loc.mrna_maxcov) continue;
      if (printfd==0) fprintf(f,"%s",loc.mrnas[i]->getID());
          else fprintf(f,",%s",loc.mrnas[i]->getID());
      printfd++;
      }
   const char* rdescr=findDescr(loc.mrna_maxcov);
   if (rdescr==NULL)  fprintf(f,"\t\n");
                 else fprintf(f,"\t%s\n",rdescr);
   }
}

void printXQ1(FILE* f, int qidx, GList<GLocus>& qloci) {
  int printfd=0;
  //print
  for (int i=0;i<qloci.Count();i++) {
     if (qloci[i]->qfidx!=qidx) continue;
      for (int j=0;j<qloci[i]->mrnas.Count();j++) {
        if (printfd==0) fprintf(f,"%s",qloci[i]->mrnas[j]->getID());
            else fprintf(f,",%s",qloci[i]->mrnas[j]->getID());
        printfd++;
        }
      }
  if (printfd==0) fprintf(f,"-");
 }

void numXLoci(GList<GXLocus>& xloci, int& last_id) {
  for (int l=0;l<xloci.Count();l++) {
    if (xloci[l]->qloci.Count()==0) continue; //we never print ref-only xloci
    last_id++;
    xloci[l]->id=last_id;
    }
}


class GProtCl {
 public:
   GList<GXConsensus> protcl;
   GProtCl(GXConsensus* c=NULL):protcl(true,false,false) {
    if (c!=NULL)
       protcl.Add(c);
    }
   bool add_Pcons(GXConsensus* c) {
    if (c==NULL || c->aalen==0) return false;
    if (protcl.Count()==0) {
        protcl.Add(c);
        return true;
        }
    if (protcl[0]->aalen!=c->aalen) return false;
    if (strcmp(protcl[0]->aa,c->aa)!=0) return false;
    protcl.Add(c);
    return true;
    }

   void addMerge(GProtCl& pcl, GXConsensus* pclnk) {
    for (int i=0;i<pcl.protcl.Count();i++) {
      if (pcl.protcl[i]!=pclnk) {
          protcl.Add(pcl.protcl[i]);
          }
      }
    }

   int aalen() {
    if (protcl.Count()==0) return 0;
    return protcl[0]->aalen;
   }
   bool operator==(GProtCl& cl) {
    return this==&cl;
    }
   bool operator<(GProtCl& cl) {
    return (this<&cl);
    }
};

class GTssCl:public GSeg { //experiment cluster of ref loci (isoforms)
 public:
   uint fstart; //lowest coordinate of the first exon
   uint fend; //highest coordinate of the first exon
   GList<GXConsensus> tsscl;
   GTssCl(GXConsensus* c=NULL):tsscl(true,false,false) {
     start=0;
     end=0;
     fstart=0;
     fend=0;
     if (c!=NULL) addFirst(c);
     }

   void addFirst(GXConsensus* c) {
     tsscl.Add(c);
     start=c->start;
     end=c->end;
     GffExon* fexon=(c->tcons->strand=='-') ? c->tcons->exons.Last() :
                                             c->tcons->exons.First();
     fstart=fexon->start;
     fend=fexon->end;
     }
   bool add_Xcons(GXConsensus* c) {
     if (tsscl.Count()==0) {
            addFirst(c);
            return true;
            }
     //check if it can be added to existing xconsensi
     uint nfend=0;
     uint nfstart=0;
     /*
     if (tsscl.Get(0)->tcons->getGeneID()!=NULL &&
             c->tcons->getGeneID()!=NULL &&
            strcmp(tsscl.Get(0)->tcons->getGeneID(), c->tcons->getGeneID()))
        //don't tss cluster if they don't have the same GeneID (?)
        //CHECKME: we might not want this if input files are not from Cufflinks
        //       and they could simply lack proper GeneID
          return false;
      */
     if (c->tcons->strand=='-') {
        //no, the first exons don't have to overlap
        //if (!c->tcons->exons.Last()->overlap(fstart,fend)) return false;
        nfstart=c->tcons->exons.Last()->start;
        nfend=c->tcons->exons.Last()->end;
        //proximity check for the transcript start:
        if (nfend>fend+tssDist || fend>nfend+tssDist)
                return false;
        }
      else {
        //if (!c->tcons->exons.First()->overlap(fstart,fend)) return false;
        nfstart=c->tcons->exons.First()->start;
        nfend=c->tcons->exons.First()->end;
        if (nfstart>fstart+tssDist || fstart>nfstart+tssDist)
            return false;
        }
     // -- if we are here, we can add to tss cluster

     tsscl.Add(c);
     if (fstart>nfstart) fstart=nfstart;
     if (fend<nfend) fend=nfend;
     if (start>c->start) start=c->start;
     if (end<c->end) end=c->end;
     return true;
     }

   void addMerge(GTssCl& cl, GXConsensus* clnk) {
     for (int i=0;i<cl.tsscl.Count();i++) {
         if (cl.tsscl[i]==clnk) continue;
         tsscl.Add(cl.tsscl[i]);
         }
     if (fstart>cl.fstart) fstart=cl.fstart;
     if (fend<cl.fend) fend=cl.fend;
     if (start>cl.start) start=cl.start;
     if (end<cl.end) end=cl.end;
     }
};
/*
class IntArray { //two dimensional int array
    int* mem;
    int xsize;
    int ysize;
  public:
   IntArray(int xlen, int ylen) {
     xsize=xlen;
     ysize=ylen;
     GMALLOC(mem, xsize*ysize*sizeof(int));
     }
   ~IntArray() {
     GFREE(mem);
     }
   int& data(int x, int y) {
    return mem[y*xsize+x];
    }
};

int aa_diff(GXConsensus* c1, GXConsensus* c2) {
 int diflen=abs(c1->aalen-c2->aalen);
 if (diflen>=6) return diflen;
 //obvious case: same CDS
 if (diflen==0 && strcmp(c1->aa, c2->aa)==0) return 0;
 //simple edit distance calculation
 IntArray dist(c1->aalen+1, c2->aalen+1);
 for (int i=0;i<=c1->aalen;i++) {
     dist.data(i,0) = i;
     }
 for (int j = 0; j <= c2->aalen; j++) {
     dist.data(0,j) = j;
     }
 for (int i = 1; i <= c1->aalen; i++)
     for (int j = 1; j <= c2->aalen; j++) {
         dist.data(i,j) = GMIN3( dist.data(i-1,j)+1,
             dist.data(i,j-1)+1,
                 dist.data(i-1,j-1)+((c1->aa[i-1] == c2->aa[j-1]) ? 0 : 1) );
         }
 int r=dist.data(c1->aalen,c2->aalen);
 return r;
}
*/
void printConsGTF(FILE* fc, GXConsensus* xc, int xlocnum) {
 GStr t_id, g_id, xloc;
 GStr gene_name(xc->tcons->getGeneName());//always preserve original gene_name if present
 if (gffAnnotate) {
	 t_id=xc->tcons->getID();
	 g_id=xc->tcons->getGeneID();
	 if (g_id.is_empty()) g_id.appendfmt("XLOC_%06d",xlocnum);
	 else xloc.appendfmt("XLOC_%06d",xlocnum);
 }
  else {
 	 t_id.appendfmt("%s_%08d",cprefix, xc->id);
	 g_id.appendfmt("XLOC_%06d",xlocnum);
 }
 fprintf(fc,
   "%s\t%s\ttranscript\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\"; gene_id \"%s\";"  ,
   xc->tcons->getGSeqName(),xc->tcons->getTrackName(),xc->tcons->start, xc->tcons->end, xc->tcons->strand,
     t_id.chars(), g_id.chars());
 GStr ref_gene_name;
 GStr ref_gene_id;
 if (gene_name.is_empty() && xc->ref!=NULL && xc->refcode>0 && classcode_rank(xc->refcode)<15) {
	//TODO: what if we want to override existing one?
    ref_gene_name=xc->ref->getGeneName();//get the gene name from this overlapping reference
    if (ref_gene_name.is_empty())
    	ref_gene_name=xc->ref->getGeneID(); //last resort: use reference gene ID (might be meaningless!)
    gene_name=ref_gene_name;
 }
 if (!gene_name.is_empty())
	 fprintf (fc, " gene_name \"%s\";", gene_name.chars());

 if (!xloc.is_empty())
	 fprintf(fc, " xloc \"%s\";",xloc.chars());

 if (gffAnnotate) {
	 //preserve ref_gene_id and ref_gene_name attributes if found, or replace them with new ones!
	 char* s=xc->tcons->getAttr("ref_gene_id", true);
	 if (s) fprintf(fc, " ref_gene_id \"%s\";", s); //TODO: what if we want to override existing one?
	 else if (haveRefs) {
		   if (ref_gene_id.is_empty() && xc->ref!=NULL && xc->refcode>0 && classcode_rank(xc->refcode)<CLASSCODE_OVL_RANK) {
		     ref_gene_id=xc->ref->getGeneID();//get the gene name from this overlapping reference
		   }
	      if (!ref_gene_id.is_empty() && ref_gene_id!=xc->tcons->getGeneID())
	        fprintf(fc, " ref_gene_id \"%s\";",ref_gene_id.chars());
	 }
	 s=xc->tcons->getAttr("ref_gene_name", true);
	 if (s) fprintf(fc, " ref_gene_name \"%s\";", s); //TODO: unless we want to override it !
	 else if (haveRefs) {
		   if (ref_gene_name.is_empty() && xc->ref!=NULL && xc->refcode>0 && classcode_rank(xc->refcode)<CLASSCODE_OVL_RANK) {
		     ref_gene_name=xc->ref->getGeneName();//get the gene name from this overlapping reference
		   }
	      if (!ref_gene_name.is_empty() && ref_gene_name!=gene_name)
	        fprintf(fc, " ref_gene_name \"%s\";",ref_gene_name.chars());
	 }
 }
 if (gffAnnotate) {
	      if (xc->contained) {
	        fprintf(fc, " contained_in \"%s\";", xc->contained->tcons->getID());
	        }
 }
 else {
     fprintf(fc, " oId \"%s\";",xc->tcons->getID());
     if (xc->contained) {
       fprintf(fc, " contained_in \"%s_%08d\";", cprefix, xc->contained->id);
       }
 }
 if (haveRefs) {
    if (xc->ref)
        fprintf(fc, " cmp_ref \"%s\";",xc->ref->getID());
    fprintf(fc, " class_code \"%c\";",xc->refcode ? xc->refcode : '.');
    if (xc->ref) {
      GStr ref_gname(xc->ref->getGeneName());
      if (ref_gname.is_empty()) ref_gname=xc->ref->getGeneID();
      if (!ref_gname.is_empty() && ref_gname!=gene_name)
        fprintf(fc, " cmp_ref_gene \"%s\";",ref_gname.chars());
    }
 }
 else { //if no reference, preserve existing class annotation if any
	char* s=xc->tcons->getAttr("cmp_ref", true);
	if (s) fprintf(fc, " cmp_ref \"%s\";",s);
	s=xc->tcons->getAttr("class_code", true);
	if (s) fprintf(fc, " class_code \"%s\";", s);
	s=xc->tcons->getAttr("cmp_ref_gene", true);
	if (s) fprintf(fc, " cmp_ref_gene \"%s\";",s);
 }
 if (xc->tss_id>0) fprintf(fc, " tss_id \"TSS%d\";",xc->tss_id);
 if (xc->p_id>0) fprintf(fc, " p_id \"P%d\";",xc->p_id);
 if (numQryFiles>1) {
	 fprintf(fc, " num_samples \"%d\";",xc->qcount);
 }
 fprintf(fc,"\n");
 //now print exons
 for (int i=0;i<xc->tcons->exons.Count();i++) {
   fprintf(fc,
     "%s\t%s\texon\t%d\t%d\t.\t%c\t.\ttranscript_id \"%s\"; gene_id \"%s\"; exon_number \"%d\";\n",
     xc->tcons->getGSeqName(),xc->tcons->getTrackName(),xc->tcons->exons[i]->start, xc->tcons->exons[i]->end, xc->tcons->strand,
	 t_id.chars(), g_id.chars(), i+1);
       //xlocnum, cprefix, xc->id, i+1);
   }
}

void tssCluster(GXLocus& xloc)
{
    GList<GTssCl> xpcls(true,true,false);
    for (int i=0;i<xloc.tcons.Count();i++)
    {
        GXConsensus* c=xloc.tcons[i];
        //if (c->tcons->exons.Count()<2) continue;  //should we skip single-exon transcripts ??
        GArray<int> mrgloci(true);
        int lfound=0;
        for (int l=0;l<xpcls.Count();l++)
        {
            if (xpcls[l]->end<c->tcons->exons.First()->start) continue;
            if (xpcls[l]->start>c->tcons->exons.Last()->end) break;
            if (xpcls[l]->add_Xcons(c))
            {
                lfound++;
                mrgloci.Add(l);

            }

        } // for each xpcluster
        if (lfound==0)
        {
            //create a xpcl with only this xconsensus
            xpcls.Add(new GTssCl(c));

        }
        else if (lfound>1)
        {
            for (int l=1;l<lfound;l++)
            {
                int mlidx=mrgloci[l]-l+1;
                xpcls[mrgloci[0]]->addMerge(*xpcls[mlidx], c);
                xpcls.Delete(mlidx);
            }
        }

    }//for each xconsensus in this xlocus
    for (int l=0;l<xpcls.Count();l++)
    {
        //if (xpcls[l]->tsscl.Count()<2) continue;
        tsscl_num++;
        for (int i=0;i<xpcls[l]->tsscl.Count();i++)
            xpcls[l]->tsscl[i]->tss_id=tsscl_num;
        //processTssCl(xcds_num, xpcls[l], faseq);
    }
}

void printXLoci(FILE* f, FILE* fc, int qcount, GList<GXLocus>& xloci, /* GFaSeqGet *faseq, */ FILE* fr=NULL) {
	for (int l=0;l<xloci.Count();l++) {
		if (xloci[l]->qloci.Count()==0) continue;
		GXLocus& xloc=*(xloci[l]);
		xloc.checkContainment(keepAltTSS, allowIntronSticking);
		tssCluster(xloc);//cluster and assign tss_id and cds_id to each xconsensus in xloc
		//protCluster(xloc,faseq);
		for (int c=0;c<xloc.tcons.Count();c++) {
			if (discardContained && xloc.tcons[c]->contained!=NULL) {
				if (!(keepRefMatching && xloc.tcons[c]->refcode=='='))
				{
					if (fr) printConsGTF(fr, xloc.tcons[c], xloc.id);
					continue;
				}
			}
			printConsGTF(fc,xloc.tcons[c],xloc.id);
			++outConsCount;
		}
		fprintf(f,"XLOC_%06d\t%s[%c]%d-%d\t", xloc.id,
				xloc.qloci[0]->mrna_maxcov->getGSeqName(),
				xloc.strand, xloc.start,xloc.end);
		//now print all transcripts in this locus, comma delimited
		//first, ref loci, if any
		int printfd=0;
		if (xloc.rloci.Count()>0) {
			for (int i=0;i<xloc.rloci.Count();i++) {
				for (int j=0;j<xloc.rloci[i]->mrnas.Count();j++) {
					if (printfd==0) fprintf(f,"%s|%s",getGeneID(xloc.rloci[i]->mrnas[j]),
							xloc.rloci[i]->mrnas[j]->getID());
					else fprintf(f,",%s|%s",getGeneID(xloc.rloci[i]->mrnas[j]),
							xloc.rloci[i]->mrnas[j]->getID());
					printfd++;
				}
			}
		}
		else {
			fprintf(f,"-");
		}
		//second, all the query transcripts
		for (int qi=0;qi<qcount;qi++) {
			fprintf(f,"\t");
			printXQ1(f,qi,xloc.qloci);
		}
		fprintf(f,"\n");
	}
}

void writeIntron(FILE* f, char strand, GFaSeqGet* faseq, GSeg& iseg,
                GList<GffObj>& mrnas, bool wrong=false) {
//find a ref mrna having this intron
  GffObj* rm=NULL;
  for (int i=0;i<mrnas.Count();i++) {
   GffObj* m=mrnas[i];
   if (m->start>iseg.end) break;
   if (m->end<iseg.start) continue;
   //intron coords overlaps mrna region
   for (int j=1;j<m->exons.Count();j++) {
      if (iseg.start==m->exons[j-1]->end+1 &&
            iseg.end==m->exons[j]->start-1) { rm=m; break; } //match found
      }//for each intron
   if (rm!=NULL) break;
   } //for each ref mrna in this locus
 if (rm==NULL) GError("Error: couldn't find ref mrna for intron %d-%d! (BUG)\n",
                         iseg.start,iseg.end);
 int ilen=iseg.end-iseg.start+1;
 fprintf(f,"%s\t%s\tintron\t%d\t%d\t.\t%c\t.\t",
            rm->getGSeqName(),rm->getTrackName(),iseg.start,iseg.end,strand);
 if (faseq!=NULL) {
   const char* gseq=faseq->subseq(iseg.start, ilen);
   char* cseq=Gstrdup(gseq, gseq+ilen-1);
   if (strand=='-') reverseComplement(cseq, ilen);
   fprintf(f,"spl \"%c%c..%c%c\"; ", toupper(cseq[0]),toupper(cseq[1]),
           toupper(cseq[ilen-2]),toupper(cseq[ilen-1]));
   GFREE(cseq);
   }
  fprintf(f,"transcript_id \"%s\";", rm->getID());
  if (wrong) fprintf(f," noOvl=1;");
  fprintf(f,"\n");

}

void reportMIntrons(FILE* fm, FILE* fn, FILE* fq, char strand,
            GList<GSuperLocus>& cmpdata) {
  if (fm==NULL) return;
  for (int l=0;l<cmpdata.Count();l++) {
    GSuperLocus *sl=cmpdata[l];
    //cache the whole locus sequence if possible
    //write these introns and their splice sites into the file
    for (int i=0;i<sl->i_missed.Count();i++)
      writeIntron(fm, strand, NULL, sl->i_missed[i], sl->rmrnas);
    if (fn!=NULL) {
        for (int i=0;i<sl->i_notp.Count();i++)
          writeIntron(fn, strand, NULL, sl->i_notp[i], sl->rmrnas);
        }
    if (fq!=NULL) {
     for (int i=0;i<sl->i_qwrong.Count();i++) {
       writeIntron(fq, strand, NULL, sl->i_qwrong[i], sl->qmrnas, true);
       }
     for (int i=0;i<sl->i_qnotp.Count();i++) {
       writeIntron(fq, strand, NULL, sl->i_qnotp[i], sl->qmrnas);
       }
     }
    }
}

void writeNIntron(FILE* f, char strand, GFaSeqGet* faseq, GSeg& iseg,
                GList<GffObj>& mrnas) {
//find all the mrnas having this intron
  GVec<GffObj*> rms;
  for (int i=0;i<mrnas.Count();i++) {
   GffObj* m=mrnas[i];
   if (m->start>iseg.end) break;
   if (m->end<iseg.start) continue;
   //intron coords overlaps mrna region
   for (int j=1;j<m->exons.Count();j++) {
      if (iseg.start==m->exons[j-1]->end+1 &&
            iseg.end==m->exons[j]->start-1) {
              rms.Add(m);
            } //match found
      }//for each intron
   } //for each ref mrna in this locus
 if (rms.Count()==0) GError("Error: couldn't find transcripts for intron %d-%d! (BUG)\n",
                         iseg.start,iseg.end);
 int ilen=iseg.end-iseg.start+1;
 fprintf(f,"%s\t%d\t%d\t%c\t",
            rms[0]->getGSeqName(),iseg.start,iseg.end,strand);

 if (faseq!=NULL) { //print splice sites!
   const char* gseq=faseq->subseq(iseg.start, ilen);
   char* cseq=Gstrdup(gseq, gseq+ilen-1);
   if (strand=='-') reverseComplement(cseq, ilen);
   fprintf(f,"%c%c..%c%c\t", toupper(cseq[0]),toupper(cseq[1]),
           toupper(cseq[ilen-2]),toupper(cseq[ilen-1]));
   GFREE(cseq);
 }
 for (int i=0;i<rms.Count();i++) {
    if (i) fprintf(f, ",%s", rms[i]->getID());
    else fprintf(f, "%s", rms[i]->getID());
 }
 fprintf(f,"\n");
}

void reportNIntrons(FILE* fn, GFaSeqGet* fq, char strand, GList<GSuperLocus>& cmpdata) {
	if (fn==NULL) return;
	for (int l=0;l<cmpdata.Count();l++) {
		GSuperLocus *sl=cmpdata[l];
		for (int i=0;i<sl->i_qnotp.Count();i++)
			writeNIntron(fn, strand, fq, sl->i_qnotp[i], sl->qmrnas);
	}
}

void processLoci(GSeqData& seqdata, GSeqData* refdata, int qfidx) {
    //GList<GSeqLoci>& glstloci, GList<GSeqCmpRegs>& cmpdata)

  if (refdata!=NULL) {
     //if (gtf_tracking_verbose) GMessage(" ..comparing to reference loci..\n") ;
     compareLoci2R(seqdata.loci_f, seqdata.gstats_f, refdata->loci_f, qfidx);
     compareLoci2R(seqdata.loci_r, seqdata.gstats_r, refdata->loci_r, qfidx);
     // -- report

     if (f_mintr!=NULL) {
       //GMessage(" ..reporting missed ref introns..\n");
       //reportIntrons(f_mintr, f_nintr, f_qintr, faseq, '+', seqdata.gstats_f);
       //reportIntrons(f_mintr, f_nintr, f_qintr, faseq, '-', seqdata.gstats_r);
       reportMIntrons(f_mintr, NULL, NULL, '+', seqdata.gstats_f);
       reportMIntrons(f_mintr, NULL, NULL, '-', seqdata.gstats_r);
       }
     }
  if (f_nj!=NULL) {
	  reportNIntrons(f_nj, NULL, '+', seqdata.gstats_f);
	  reportNIntrons(f_nj, NULL, '-', seqdata.gstats_r);
  }
}

void collectRLocData(GSuperLocus& stats, GLocus& loc) {
	//this is only called for ref loci not overlapping any query loci
	stats.total_rloci++;
	stats.total_rmrnas+=loc.mrnas.Count();
	stats.total_richains+=loc.ichains;
	stats.total_rmexons+=loc.mexons.Count();
	stats.total_rexons+=loc.uexons.Count();
	stats.total_rintrons+=loc.introns.Count();
	stats.m_exons+=loc.uexons.Count(); //missed ref exons
	stats.m_introns+=loc.introns.Count(); //missed ref introns
	stats.m_loci++; //missed ref loci
	for (int e=0;e<loc.mexons.Count();e++) {
		 stats.rbases_all+=loc.mexons[e].end-loc.mexons[e].start+1;
	}
}

void collectRNOvl(GSuperLocus& stats, GList<GLocus>& loci, int qfidx) { //, const char* gseqname) {
  for (int l=0;l<loci.Count();l++) {
    if (loci[l]->qlocovls[qfidx]->Count()==0) // cmpovl is across all query files so it won't account
    	                            //for refs if a previous query set cmpovl for that ref locus
      collectRLocData(stats,*loci[l]);
  }
}

void collectRData(GSuperLocus& stats, GList<GLocus>& loci) {
 for (int l=0;l<loci.Count();l++)
      collectRLocData(stats,*loci[l]);
}

//adjust stats for a list of unoverlapped (completely "wrong" or novel) qry loci
void collectQLocData(GSuperLocus& stats, GLocus& loc) {
 stats.total_qmrnas+=loc.mrnas.Count();
 stats.total_qexons+=loc.uexons.Count();
 stats.total_qmexons+=loc.mexons.Count();
 stats.total_qintrons+=loc.introns.Count();
 stats.total_qichains+=loc.ichains;
 stats.total_qloci++;
 stats.w_loci++; //add to the count of novel/wrong loci
 if (loc.ichains>0 && loc.mrnas.Count()>1)
    stats.total_qloci_alt++;
 stats.w_exons+=loc.uexons.Count();
 stats.w_introns+=loc.introns.Count();
 for (int e=0;e<loc.mexons.Count();e++) {
   stats.qbases_all+=loc.mexons[e].end-loc.mexons[e].start+1;
   }
}

void collectQData(GSuperLocus& stats, GList<GLocus>& loci) {
 for (int l=0;l<loci.Count();l++) {
     //this is called when no refdata is given, so all these loci are nloci
     //nloci.Add(loci[l]);
     collectQLocData(stats,*loci[l]);
     }
}

void collectQNOvl(GSuperLocus& stats, GList<GLocus>& loci, GList<GLocus>& nloci) {
  for (int l=0;l<loci.Count();l++) {
    if (loci[l]->cmpovl.Count()==0) {//locus with no ref loci overlaps
      nloci.Add(loci[l]);
      if (!reduceQrys) collectQLocData(stats,*loci[l]);
      }
  }
}

void collectQU(GSuperLocus& stats, GList<GLocus>& nloci) {
  for (int l=0;l<nloci.Count();l++) {
    //stats.w_loci++; //novel/wrong loci
    collectQLocData(stats, *nloci[l]);
    }
}

void printLocus(FILE* f, GLocus& loc, const char* gseqname) {
  fprintf(f, "## Locus %s:%d-%d\n",gseqname, loc.start, loc.end);
  for (int m=0;m<loc.mrnas.Count();m++) {
    loc.mrnas[m]->printGtf(f);
    }
}

void collectCmpData(GSuperLocus& stats, GList<GSuperLocus>& cmpdata) { //, const char* gseqname) {
 for (int c=0;c<cmpdata.Count();c++) {
   stats.addStats(*cmpdata[c]);
   /*
   if (f_nloci!=NULL && cmpdata[c]->locusTP==0 && cmpdata[c]->rloci.Count()>0) {
      fprintf(f_nloci, "# Superlocus %s:%d-%d\n",gseqname, cmpdata[c]->start, cmpdata[c]->end);
      for (int l=0;l<cmpdata[c]->rloci.Count();l++) {
         printLocus(f_nloci,*cmpdata[c]->rloci[l], gseqname);
         }
      }
   */
   }
}

void printLociQ(FILE* f, GList<GLocus>& loci, char c=' ') {
  for (int i=0;i<loci.Count();i++) {
     GLocus& loc=*(loci[i]);
     for (int j=0;j<loc.mrnas.Count();j++) {
        fprintf(f, "%s\t%c\t%c\t%s\n", loc.mrnas[j]->getID(), loc.mrnas[j]->strand, c, loc.mrnas[j]->getGSeqName());
     }
  }
}

void printLocQ(FILE *f, GLocus& loc) {
    for (int j=0;j<loc.mrnas.Count();j++) {
       fprintf(f, "%s\n", loc.mrnas[j]->getID());
    }
}


void collectStats(GSuperLocus& stats, GSeqData* seqdata, GSeqData* refdata, int qfidx) {
 //collect all stats for a single genomic sequence into stats
 if (seqdata==NULL) { //for contig/chromosome with no qry data
   if (reduceRefs || refdata==NULL) return;
   collectRData(stats, refdata->loci_f);
   collectRData(stats, refdata->loci_r);
   return;
 }

 if (refdata==NULL) {//reference data missing on this contig
   for (int l=0;l<seqdata->loci_f.Count();l++) {
      seqdata->nloci_f.Add(seqdata->loci_f[l]);
   }
   for (int l=0;l<seqdata->loci_r.Count();l++) {
      seqdata->nloci_r.Add(seqdata->loci_r[l]);
   }
   if (reduceQrys) return;

   collectQData(stats, seqdata->loci_f);
   collectQData(stats, seqdata->loci_r);
   collectQU(stats, seqdata->nloci_u);
   return;
 }
 //collect data for overlapping superloci (already in seqdata->gstats_f/_r)
 collectCmpData(stats, seqdata->gstats_f);
 collectCmpData(stats, seqdata->gstats_r);
 //for non-overlapping qry loci, add them as false positives FP

 collectQNOvl(stats, seqdata->loci_f , seqdata->nloci_f);
 collectQNOvl(stats, seqdata->loci_r , seqdata->nloci_r);
 if (!reduceQrys) {
   collectQU(stats, seqdata->nloci_u);
 }
 if (!reduceRefs) {
   collectRNOvl(stats, refdata->loci_f, qfidx);
   collectRNOvl(stats, refdata->loci_r, qfidx);
 }
}

void reportStats(FILE* fout, const char* setname, GSuperLocus& stotal,
                          GSeqData* seqdata, GSeqData* refdata, int qfidx) {
  GSuperLocus stats;
  bool finalSummary=(seqdata==NULL && refdata==NULL);
  GSuperLocus *ps=(finalSummary ? &stotal : &stats );
  if (!finalSummary) { //collecting contig stats
    //gather statistics for all loci/superloci here
    collectStats(stats, seqdata, refdata, qfidx);
    if (f_qdisc!=NULL) {
    	printLociQ(f_qdisc, seqdata->nloci_f, 'F');
    	printLociQ(f_qdisc, seqdata->nloci_r, 'R');
    	printLociQ(f_qdisc, seqdata->nloci_u, 'U');
    }
    stotal.addStats(stats);
    if (!perContigStats) return;
  }

  ps->calcF();
  if (seqdata!=NULL) fprintf(fout, "#> Genomic sequence: %s \n", setname);
                else fprintf(fout, "\n#= Summary for dataset: %s \n", setname);

  fprintf(fout,   "#     Query mRNAs : %7d in %7d loci  (%d multi-exon transcripts)\n",
          ps->total_qmrnas, ps->total_qloci, ps->total_qichains);
  fprintf(fout, "#            (%d multi-transcript loci, ~%.1f transcripts per locus)\n",
          ps->total_qloci_alt, ((double)ps->total_qmrnas/ps->total_qloci));

  if (haveRefs) {
    fprintf(fout, "# Reference mRNAs : %7d in %7d loci  (%d multi-exon)\n",
            ps->total_rmrnas, ps->total_rloci, ps->total_richains);
    if (ps->baseTP+ps->baseFP==0 || ps->baseTP+ps->baseFN==0) return;
    fprintf(fout, "# Super-loci w/ reference transcripts:  %7d\n",ps->total_superloci);

    /*if (seqdata!=NULL) {
      fprintf(fout, "          ( %d/%d on forward/reverse strand)\n",
             seqdata->gstats_f.Count(),seqdata->gstats_r.Count());
       }*/
    fprintf(fout, "#-----------------| Sensitivity | Precision  |\n");
    double sp=(100.0*(double)ps->baseTP)/(ps->baseTP+ps->baseFP);
    double sn=(100.0*(double)ps->baseTP)/(ps->baseTP+ps->baseFN);
    fprintf(fout, "        Base level:   %5.1f     |   %5.1f    |\n",sn, sp);
    sp=(100.0*(double)ps->exonQTP)/ps->total_qexons;
    sn=(100.0*(double)ps->exonTP)/ps->total_rexons;
    //DEBUG only:
    //fprintf(fout, "======> Exon stats: %d total_qexons, %d total_rexons, %d exonTP, %d exonFP, %d exonFN\n",
    //   ps->total_rexons, ps->total_qexons, ps->exonTP, ps->exonFP, ps->exonFN);
    fprintf(fout, "        Exon level:   %5.1f     |   %5.1f    |\n",sn, sp);
    if (ps->total_rintrons>0) {
      //intron level
      sp=(100.0*(double)ps->intronTP)/(ps->intronTP+ps->intronFP);
      sn=(100.0*(double)ps->intronTP)/(ps->intronTP+ps->intronFN);
      fprintf(fout, "      Intron level:   %5.1f     |   %5.1f    |\n",sn, sp);
      //intron chains:
      sp=(100.0*(double)ps->ichainTP)/ps->total_qichains;
      sn=(100.0*(double)ps->ichainTP)/ps->total_richains;
      fprintf(fout, "Intron chain level:   %5.1f     |   %5.1f    |\n",sn, sp);
    }
    sp=(100.0*(double)ps->mrnaTP)/ps->total_qmrnas;
    sn=(100.0*(double)ps->mrnaTP)/ps->total_rmrnas;
    fprintf(fout, "  Transcript level:   %5.1f     |   %5.1f    |\n",sn, sp);
    sp=(100.0*(double)ps->locusQTP)/ps->total_qloci;
    sn=(100.0*(double)ps->locusTP)/ps->total_rloci;  //(ps->locusTP+ps->locusFN);
    fprintf(fout, "       Locus level:   %5.1f     |   %5.1f    |\n",sn, sp);
    //fprintf(fout, "                   (locus TP=%d, total ref loci=%d)\n",ps->locusTP, ps->total_rloci);
    fprintf(fout,"\n     Matching intron chains: %7d\n",ps->ichainTP);
    fprintf(fout,  "       Matching transcripts: %7d\n",ps->mrnaTP);
    fprintf(fout,  "              Matching loci: %7d\n",ps->locusTP);
    fprintf(fout, "\n");
    sn=(100.0*(double)ps->m_exons)/(ps->total_rexons);
    fprintf(fout, "          Missed exons: %7d/%d\t(%5.1f%%)\n",ps->m_exons, ps->total_rexons, sn);
    sn=(100.0*(double)ps->w_exons)/(ps->total_qexons);
    fprintf(fout, "           Novel exons: %7d/%d\t(%5.1f%%)\n",ps->w_exons, ps->total_qexons,sn);
    if (ps->total_rintrons>0) {
    sn=(100.0*(double)ps->m_introns)/(ps->total_rintrons);
    fprintf(fout, "        Missed introns: %7d/%d\t(%5.1f%%)\n",ps->m_introns, ps->total_rintrons, sn);
    }
    if (ps->total_qintrons>0) {
    sn=(100.0*(double)ps->w_introns)/(ps->total_qintrons);
    fprintf(fout, "         Novel introns: %7d/%d\t(%5.1f%%)\n",ps->w_introns, ps->total_qintrons,sn);
    }
    if (ps->total_rloci>0) {
    sn=(100.0*(double)ps->m_loci)/(ps->total_rloci);
    fprintf(fout, "           Missed loci: %7d/%d\t(%5.1f%%)\n",ps->m_loci, ps->total_rloci, sn);
    }
    if (ps->total_qloci>0) {
    sn=(100.0*(double)ps->w_loci)/(ps->total_qloci);
    fprintf(fout, "            Novel loci: %7d/%d\t(%5.1f%%)\n",ps->w_loci, ps->total_qloci,sn);
    }

  }
}

int inbuf_len=1024; //starting inbuf capacity
char* inbuf=NULL; // incoming buffer for sequence lines.

void loadRefDescr(const char* fname) {
  if (inbuf==NULL)  { GMALLOC(inbuf, inbuf_len); }
  FILE *f=fopen(fname, "rb");
  if (f==NULL) GError("Error opening exon file: %s\n",fname);
  char* line;
  int llen=0;
  off_t fpos;
  while ((line=fgetline(inbuf, inbuf_len, f, &fpos, &llen))!=NULL) {
   if (strlen(line)<=2) continue;
   int idlen=strcspn(line,"\t ");
   char* p=line+idlen;
   if (idlen<llen && idlen>0) {
     *p=0;
      p++;
      refdescr.Add(line, new GStr(p));
      }
  }
}

GSeqTrack* findGSeqTrack(int gsid) {
  GSeqTrack f(numQryFiles, gsid);
  int fidx=-1;
  if (gseqtracks.Found(&f,fidx))
     return gseqtracks[fidx];
  fidx=gseqtracks.Add(new GSeqTrack(numQryFiles, gsid));
  return gseqtracks[fidx];
}

GffObj* findRefMatch(GffObj& m, GLocus& rloc, int& ovlen) {
 ovlen=0;
 CTData* mdata=((CTData*)m.uptr);
 if (mdata->eqref!=NULL && ((CTData*)(mdata->eqref->uptr))->locus==&rloc) {
	  if (mdata->eqref!=mdata->ovls.First()->mrna) {
		  GMessage("Warning: previously matching ref (%s) not the best overlap for %s (which is %s)\n",
				  mdata->eqref->getID(), m.getID(), mdata->ovls.First()->mrna->getID());
		  mdata->eqref=mdata->ovls.First()->mrna;
	  }

      return mdata->eqref;
 }
 //if (rloc==NULL|| m==NULL) return NULL;
 GffObj* ret=NULL;
 for (int r=0;r<rloc.mrnas.Count();r++) {
    int olen=0;
    char eqcode=0;
	if ((eqcode=transcriptMatch(m, *(rloc.mrnas[r]),olen))>0) {
	  /*
	  if (ovlen<olen) {
		  ovlen=olen;
		  ret=rloc.mrnas[r]; //keep the longest matching ref
			 //but this is unnecessary, there can be only one matching ref
			 // because duplicate refs were discarded
		  }
	  */
	  //for class code output, '~' should be shown as '=' unless strict matching was requested!
	  if (eqcode=='~' && !stricterMatching) { eqcode='='; olen--; }
	  mdata->addOvl(eqcode,rloc.mrnas[r], olen);
      //this must be called only for the head of an equivalency chain
      CTData* rdata=(CTData*)rloc.mrnas[r]->uptr;
      rdata->addOvl(eqcode,&m,olen);
      //if (rdata->eqnext==NULL) rdata->eqnext=&m;
      }
    }
 if (mdata->ovls.Count()>0 && mdata->ovls.First()->code=='=') {
   ret=mdata->ovls.First()->mrna;
   ovlen=mdata->ovls.First()->ovlen;
 }
 //
 if (ret!=NULL)
   mdata->eqref=ret;
 return ret;
}


GXConsensus* addXCons(GXLocus* xloc, GffObj* ref, char ovlcode, GffObj* tcons, CEqList* ts) {
 GXConsensus* c=new GXConsensus(tcons, ts, ref, ovlcode);
 xloc->tcons.Add(c);
 return c;
}

void findTMatches(GTrackLocus& loctrack, int qcount) {
	//perform an all vs. all ichain-match for all transcripts across all loctrack[i]->qloci
	//link best "match" for each qry transcript across qry datasets
	for (int q=0;q<qcount-1;q++) { //for each qry dataset
		if (loctrack[q]==NULL) continue;
		for (int qi=0;qi<loctrack[q]->Count();qi++) { // for each transcript in q dataset
			GffObj* qi_t=loctrack[q]->Get(qi);
			CTData* qi_d=(CTData*)qi_t->uptr;
			//if (qi_d->eqlist!=NULL && qi_t->exons.Count()>1)
			//  continue; //this is part of an EQ chain already
			GArray<CEqMatch> eqmatches(true);
			for (int n=q+1;n<qcount;n++) { // for every next/successor dataset
				if (loctrack[n]==NULL) continue;
				for (int ni=0;ni<loctrack[n]->Count();ni++) { //for every transcript in this next dataset
					GffObj* ni_t=loctrack[n]->Get(ni);
					CTData* ni_d=(CTData*)ni_t->uptr;
					if (ni_d->eqlist!=NULL && ni_d->eqlist==qi_d->eqlist) continue;
					int ovlen=0;
					if (transcriptMatch(*qi_t, *ni_t, ovlen)>0) {
						CEqMatch m(ni_t, tMatchScore(ovlen, ni_t, qi_t));
						eqmatches.Add(&m);
					}
				}
	            //select the best match in this next dataset
                if (eqmatches.Count()>0) {
                	CEqMatch& m=eqmatches.Last();
                	qi_d->joinEqList(m.t);
                }
			} // for each successor dataset
		} //for each transcript in qry dataset
	} //for each qry dataset
}


int cmpTData_qset(const pointer* p1, const pointer* p2) {
 CTData* d1=(CTData*)(((GffObj*)p1)->uptr);
 CTData* d2=(CTData*)(((GffObj*)p2)->uptr);
 return (d1->qset - d2->qset);
 }

void printITrack(FILE* ft, GList<GffObj>& mrnas, int qcount, int& cnum) {
	for (int i=0;i<mrnas.Count();i++) {
		GffObj& qt=*(mrnas[i]);
		CTData* qtdata=(CTData*)qt.uptr;
		int qfidx=qtdata->qset;
		char ovlcode=qtdata->classcode;
		CEqList* eqchain=qtdata->eqlist;
		GffObj* ref=NULL; //related ref -- it doesn't have to be fully matching
		//GffObj* eqref=NULL; //fully ichain-matching ref
		//GffObj* tcons=NULL; //"consensus" (largest) transcript for a clique
		//int tmaxcov=0;
		//eqchain.Add(&qt);
		//eqref=qtdata->eqref;
		if (qtdata->ovls.Count()>0 && qtdata->ovls[0]->mrna!=NULL) {
			//if it has ovlcode with a ref
			ref=qtdata->ovls[0]->mrna;
			//consistency check: qtdata->ovls[0]->code==ovlcode
			// -- let tcons be a transfrag, not a ref transcript
			//tcons=eqref;
			//if (tcons!=NULL) tmaxcov=tcons->covlen;
		}
		GffObj* tcons=mrnas[i];
		int tmaxcov=tcons->covlen;
		ovlcode=qtdata->getBestCode();

		if (qtdata->eqhead) {//head of a equivalency chain
			//check if all transcripts in this chain have the same ovlcode
			//GffObj* tcons_bycode=tcons;
			//bool ovlcode_change=false;
			for (int k=0;k<qtdata->eqlist->Count();k++) {
				GffObj* m=qtdata->eqlist->Get(k);
				if (m->covlen>tmaxcov) {
					tmaxcov=m->covlen;
					tcons=m;
					ovlcode=((CTData*)m->uptr)->getBestCode();
					if (((CTData*)m->uptr)->ovls.Count()>0)
						ref=((CTData*)m->uptr)->ovls[0]->mrna;
				}
			}
			//if (ovlcode_change && tcons!=tcons_bycode)
			//	tcons=tcons_bycode;
		}//chain check
		//if (ovlcode=='p') ref=NULL; //ignore polymerase runs?
		if (ovlcode==0 || ovlcode=='-' || ovlcode=='.') {
			ovlcode = (ref==NULL) ? 'u' : '.'; //should never be '.' here
		}
		//-- print columns 1 and 2 as LOCUS_ID and TCONS_ID
		//bool chainHead=(qtdata->eqnext!=NULL && ((qtdata->eqdata & EQHEAD_TAG)!=0));
		bool chainHead=qtdata->eqhead;
		//bool noChain=((qtdata->eqdata & EQCHAIN_TAGMASK)==0);
		bool noChain=(eqchain==NULL);
		GXConsensus* xtcons=NULL;
		if (chainHead || noChain) {
			cnum++;
			if (ft!=NULL) fprintf(ft,"%s_%08d\t",cprefix,cnum);
			GXLocus* xloc=qtdata->locus->xlocus;
			if (xloc!=NULL) {
				if (ft!=NULL) fprintf(ft, "XLOC_%06d\t",xloc->id);
				if (tcons->exons.Count()>1) {
					//! only multi-exon mRNAs are counted for multi-transcript xloci !
					xloc->num_mtcons++;
					if (xloc->num_mtcons==2)
						total_xloci_alt++;
				}
			}
			else {
				//should NEVER happen!
				int fidx=qtdata->qset;
				GError("Error: no XLocus created for transcript %s (file %s) [%d, %d], on %s%c:%d-%d\n", qt.getID(),
						qryfiles[qtdata->locus->qfidx]->chars(), qtdata->locus->qfidx, fidx, qt.getGSeqName(), qt.strand, qt.start, qt.end);
			}
			xtcons=addXCons(xloc, ref, ovlcode, tcons, eqchain);
		} // if chain head or uniq entry (not part of a chain)
		if (ft==NULL) continue;
		if (chainHead) {
			//this is the start of an equivalence class as a printing chain
			if (ref!=NULL) fprintf(ft,"%s|%s\t%c", getGeneID(ref),ref->getID(), ovlcode);
			else fprintf(ft,"-\t%c", ovlcode);
			GffObj* m=mrnas[i];
			CTData* mdata=(CTData*)m->uptr;

			int lastpq=-1;
			eqchain->setUnique(false);
			// sort eqchain by qry# so we can catch duplicate/redundant transfrags in the same qry set
			eqchain->setSorted((GCompareProc*) cmpTData_qset);
			xtcons->qcount=0; //should already be 0 for chain heads

			for (int k=0;k<eqchain->Count();k++) {
				m=eqchain->Get(k);
				mdata=(CTData*)m->uptr;
				if (mdata->qset==lastpq) {
					//when a qry file has duplicates/redundant transfrags
					fprintf(ft,",%s|%s|%d|%8.6f|%8.6f|%8.6f|%d", getGeneID(m), m->getID(),
							//iround(m->gscore/10),
							m->exons.Count(),
							mdata->FPKM, mdata->TPM, mdata->cov, m->covlen);
					continue;
				}
				xtcons->qcount++;
				for (int ptab=mdata->qset-lastpq;ptab>0;ptab--)
					if (ptab>1) fprintf(ft,"\t-");
					else fprintf(ft,"\t");
				lastpq = mdata->qset;
				fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%d", lastpq+1, getGeneID(m), m->getID(),
						//iround(m->gscore/10),
						m->exons.Count(),
						mdata->FPKM, mdata->TPM, mdata->cov, m->covlen);
			}
			for (int ptab=qcount-lastpq-1;ptab>0;ptab--)
				fprintf(ft,"\t-");
			fprintf(ft,"\n");
			continue;
		} //start of eq class (printing chain)

		if (eqchain!=NULL) continue; //part of a matching chain, dealt with previously

		//--------- not in an ichain-matching class, print as singleton

		if (ref!=NULL) fprintf(ft,"%s|%s\t%c",getGeneID(ref), ref->getID(), ovlcode);
		else fprintf(ft,"-\t%c",ovlcode);
		for (int ptab=qfidx;ptab>=0;ptab--)
			if (ptab>0) fprintf(ft,"\t-");
			else fprintf(ft,"\t");
		fprintf(ft,"q%d:%s|%s|%d|%8.6f|%8.6f|%8.6f|%d",qfidx+1, getGeneID(qt), qt.getID(),
				//iround(qt.gscore/10),
				qt.exons.Count(),
				qtdata->FPKM, qtdata->TPM, qtdata->cov, qt.covlen);
		for (int ptab=qcount-qfidx-1;ptab>0;ptab--)
			fprintf(ft,"\t-");
		fprintf(ft,"\n");
	} //for each transcript
}


void findTRMatch(GTrackLocus& loctrack, int qcount, GLocus& rloc) {
	//requires loctrack to be already populated with overlapping qloci by findTMatches()
	// which also found (and tagged) all matching qry transcripts
	for (int q=0;q<qcount;q++) { //for each qry dataset
		if (loctrack[q]==NULL) continue;
		for (int qi=0;qi<loctrack[q]->Count();qi++) { // for each transcript in q dataset
			//if (loctrack[q]->cl[qi]->exons.Count()<2) continue; //skip single-exon transcripts
			GffObj& qt=*(loctrack[q]->Get(qi));
			CTData* qtdata=(CTData*)qt.uptr;
			GffObj* rmatch=NULL; //== ref match for this row
			int rovlen=0;
			//if (qtdata->eqnext!=NULL && ((qtdata->eqdata & EQHEAD_TAG)!=0)) {
			if (qtdata->eqhead) {
				//EQ chain head -- transfrag equivalency list starts here
				if (qtdata->eqref==NULL) { //find rloc overlap
					if (qt.overlap(rloc.start, rloc.end)) {
						rmatch=findRefMatch(qt, rloc, rovlen);
					}
				} else rmatch=qtdata->eqref;
				if (rmatch!=NULL) {
					//assign the same refmatch to all
					//FIXME: transitivity is NOT guaranteed for single-exon transcripts!
					for (int k=0;k<qtdata->eqlist->Count();k++) {
						GffObj* m=qtdata->eqlist->Get(k);
						((CTData*)m->uptr)->addOvl('=',rmatch,rovlen);
						continue;
					}
				}
				//if (rmatch!=NULL) continue;
			} //equivalence class (chain of intron-matching)
			//if ((qtdata->eqdata & EQCHAIN_TAGMASK)!=0) continue; //part of a matching chain, dealt with previously
			//--------- qry mrna not in a '=' matching clique
			if (qtdata->eqref==NULL) { //find any rloc overlap -- class code
				if (qt.overlap(rloc.start, rloc.end)) {
					rmatch=findRefMatch(qt, rloc, rovlen);
					if (rmatch==NULL) //&& ((CTData*)qt.uptr)->ovls.Count()==0)
						{
						//not an ichain match, look for other codes
						gatherRefLocOvls(qt, rloc);
					}
				}
			}
			else rmatch=qtdata->eqref;
		} //for each qry transcript
	}//for each qry dataset
}


bool inPolyRun(char strand, GffObj& m, GList<GLocus>* rloci, int& rlocidx) {
 //we are only here if there is no actual overlap between m and any locus in rloci
  if (rloci==NULL || rloci->Count()==0) return false; // || m.exons.Count()>1
  if (strand=='-') {
        rlocidx=qsearch_loci(m.end, *rloci);
        //returns index of locus starting just ABOVE m.end
        // or -1 if last locus start <= m.end
        GLocus* rloc=NULL;
        if (rlocidx<0) return false;
        while (rlocidx<rloci->Count()) {
           rloc=rloci->Get(rlocidx);
           if (rloc->start>m.end+polyrun_range) break;
           if (rloc->start+6>m.end) return true;
           rlocidx++;
           }
  } else { // strand == '+' (or '.' ?)
        rlocidx=qsearch_loci(m.end, *rloci);
        GLocus* rloc=NULL;
        //returns index of closest locus starting ABOVE m.end
        // or -1 if last locus start <= m.end
        if (rlocidx<0) rlocidx=rloci->Count(); //this may actually start below m.end
        while ((--rlocidx)>=0) {
          rloc=rloci->Get(rlocidx);
          if (m.start>rloc->start+GFF_MAX_LOCUS) break;
          if (m.start+6>rloc->end && m.start<rloc->end+polyrun_range) return true;
        }
  }
  return false;
}

CTData* getBestOvl(GffObj& m) {
 //CTData* mdata=(CTData*)m.uptr;
 //return mdata->getBestCode();
  if ( ((CTData*)m.uptr)->ovls.Count()>0)
     return (CTData*)m.uptr;
  return NULL;
}

void reclass_XStrand(GList<GffObj>& mrnas, GList<GLocus>* rloci) {
  //checking for relationship with ref transcripts on opposite strand
  if (rloci==NULL || rloci->Count()<1) return;
  int j=0;//current rloci index
  for (int i=0;i<mrnas.Count();i++) {
     GffObj& m=*mrnas[i];
     char ovlcode=((CTData*)m.uptr)->getBestCode();
     if (ovlcode>47 && strchr("=cjeo",ovlcode)!=NULL) continue;
     GLocus* rloc=rloci->Get(j);
     if (rloc->start>m.end) continue; //check next transfrag
     while (m.start>rloc->end && j+1<rloci->Count()) {
           j++;
           rloc=rloci->Get(j);
           }
     if (rloc->start>m.end) continue; //check next transfrag
     //m overlaps rloc:
     //check if m has a fuzzy intron overlap -> 's' (shadow, mapped on the wrong strand)
     //  then if m is contained within an intron -> 'i'
     //  otherwise it's just a plain cross-strand overlap: 'x'
     int jm=0;
     do { //while rloci overlap this transfrag (m)
       rloc=rloci->Get(j+jm);
       bool is_shadow=false;
       GffObj* sovl=NULL;
       bool is_intraintron=false;
       GffObj* iovl=NULL;
       if (rloc->introns.Count()>0) {
           for (int n=0;n<rloc->introns.Count();n++) {
              GISeg& rintron=rloc->introns[n];
              if (rintron.start>m.end) break;
              if (m.start>rintron.end) continue;
              //overlap between m and intron
              if (m.end<=rintron.end && m.start>=rintron.start) {
                  is_intraintron=true;
                  if (iovl==NULL || iovl->covlen<rintron.t->covlen) iovl=rintron.t;
                  continue;
                  }
              //check if any intron of m has a fuzz-match with rintron
              for (int e=1;e<m.exons.Count();e++) {
                 GSeg mintron(m.exons[e-1]->end+1,m.exons[e]->start-1);
                 if (rintron.coordMatch(&mintron,10)) {
                    is_shadow=true;
                    if (sovl==NULL || sovl->covlen<rintron.t->covlen) sovl=rintron.t;
                    break;
                    }
                 } //for each m intron
              } //for each intron of rloc
           }//rloc has introns
       bool xcode=true;
       if (is_shadow) { ((CTData*)m.uptr)->addOvl('s', sovl); xcode=false; }
             // else
       if (ovlcode!='i' && is_intraintron) { ((CTData*)m.uptr)->addOvl('i', iovl); xcode=false; }
       if (xcode) {
               // just plain overlap, find the overlapping mrna in rloc
               GffObj* maxovl=NULL;
               int ovlen=0;
               GffObj* max_lovl=NULL; //max len ref transcript
                       // having no exon overlap but simply range overlap (interleaved exons)
               for (int ri=0;ri<rloc->mrnas.Count();ri++) {
                  if (!m.overlap(*(rloc->mrnas[ri]))) continue;
                  int o=m.exonOverlapLen(*(rloc->mrnas[ri]));
                  if (o>0) {
                     if (o>ovlen) {
                        ovlen=o;
                        maxovl=rloc->mrnas[ri];
                        }
                     }
                    else { //no exon overlap, but still overlapping (interleaved exons)
                     if (max_lovl==NULL || max_lovl->covlen<rloc->mrnas[ri]->covlen)
                         max_lovl=rloc->mrnas[ri];
                     }
                  }
               if (maxovl) ((CTData*)m.uptr)->addOvl('x',maxovl);
                 else if (max_lovl) ((CTData*)m.uptr)->addOvl('x',max_lovl);
               } //'x'
       jm++;
       } while (j+jm<rloci->Count() && rloci->Get(j+jm)->overlap(m));
     } //for each transfrag
}

void reclass_mRNAs(char strand, GList<GffObj>& mrnas, GList<GLocus>* rloci, GFaSeqGet *faseq) {
  int rlocidx=-1;
  for (int i=0;i<mrnas.Count();i++) {
    GffObj& m=*mrnas[i];
    char ovlcode=((CTData*)m.uptr)->getBestCode();
    //if (ovlcode=='u' || ovlcode=='i' || ovlcode==0) {
    if (ovlcode=='u' || ovlcode<47) {
      //check for overlaps with ref transcripts on the other strand
      if (m.exons.Count()==1 && inPolyRun(strand, m, rloci, rlocidx)) {
         ((CTData*)m.uptr)->addOvl('p', rloci->Get(rlocidx)->mrna_maxcov);
         }
      else { //check for repeat content
         if (faseq!=NULL) {
            int seqlen;
            char* seq=m.getSpliced(faseq, false, &seqlen);
            //get percentage of lowercase
            int numlc=0;
            for (int c=0;c<seqlen;c++) if (seq[c]>='a') numlc++;
            if (numlc > seqlen/2)
               ((CTData*)m.uptr)->addOvl('r', NULL, numlc);
            GFREE(seq);
            }
         }
      } //for unassigned class
  }//for each mrna

}

/* void reclassLoci(char strand, GList<GLocus>& qloci, GList<GLocus>* rloci, GFaSeqGet *faseq) {
  for (int ql=0;ql<qloci.Count();ql++) {
    reclass_mRNAs(strand, qloci[ql]->mrnas, rloci, faseq);
    //find closest upstream ref locus for this q locus
  } //for each locus
}*/

//for a single genomic sequence, all qry data and ref data is stored in gtrack
//check for all 'u' transfrags if they are repeat ('r') or polymerase run 'p' or anything else
void umrnaReclass(int qcount,  GSeqTrack& gtrack, FILE** ftr, GFaSeqGet* faseq=NULL) {
    for (int q=0;q<qcount;q++) {
        if (gtrack.qdata[q]==NULL) continue; //no transcripts in this q dataset for this genomic seq
        //reclassLoci('+', gtrack.qdata[q]->loci_f, gtrack.rloci_f, faseq);
        //reclassLoci('-', gtrack.qdata[q]->loci_r, gtrack.rloci_r, faseq);
        reclass_mRNAs('+', gtrack.qdata[q]->mrnas_f, gtrack.rloci_f, faseq);
        reclass_mRNAs('-', gtrack.qdata[q]->mrnas_r, gtrack.rloci_r, faseq);
        reclass_mRNAs('+', gtrack.qdata[q]->umrnas, gtrack.rloci_f, faseq);
        reclass_mRNAs('-', gtrack.qdata[q]->umrnas, gtrack.rloci_r, faseq);
        //and also check for special cases with cross-strand overlaps:
        reclass_XStrand(gtrack.qdata[q]->mrnas_f, gtrack.rloci_r);
        reclass_XStrand(gtrack.qdata[q]->mrnas_r, gtrack.rloci_f);
        // print all tmap data here:
        for (int i=0;i<gtrack.qdata[q]->tdata.Count();i++) {
            CTData* mdata=gtrack.qdata[q]->tdata[i];
            if (mdata->mrna==NULL) continue; //invalidated -- removed earlier
            //GLocus* rlocus=NULL;
            mdata->classcode='u';
            GffObj* ref=NULL;
            if (mdata->ovls.Count()>0) {
                mdata->classcode=mdata->ovls[0]->code;
                ref=mdata->ovls[0]->mrna;
            }
            //if (mdata->classcode<33) mdata->classcode='u';
            if (mdata->classcode<47) mdata->classcode='u'; // if 0, '-' or '.'
            if (tmapFiles) {
                char ref_match_len[2048];
                if (ref!=NULL) {
                    sprintf(ref_match_len, "%d",ref->covlen);
                    fprintf(ftr[q],"%s\t%s\t",getGeneID(ref),ref->getID());
                    //rlocus=((CTData*)(ref->uptr))->locus;
                }
                else {
                    fprintf(ftr[q],"-\t-\t");
                    strcpy(ref_match_len, "-");
                }
                //fprintf(ftr[q],"%c\t%s\t%d\t%8.6f\t%8.6f\t%d\n", ovlcode, mdata->mrna->getID(),
                //    iround(mdata->mrna->gscore/10), mdata->FPKM, mdata->cov, mdata->mrna->covlen);
                const char* mlocname = (mdata->locus!=NULL) ? mdata->locus->mrna_maxcov->getID() : mdata->mrna->getID();
                fprintf(ftr[q],"%c\t%s\t%s\t%d\t%8.6f\t%8.6f\t%8.6f\t%d\t%s\t%s\n", mdata->classcode, getGeneID(mdata->mrna), mdata->mrna->getID(),
                        //iround(mdata->mrna->gscore/10),
                		mdata->mrna->exons.Count(),
						mdata->FPKM, mdata->TPM, mdata->cov, mdata->mrna->covlen, mlocname, ref_match_len);
            }
        } //for each tdata
    } //for each qdata
}

void buildXLoci(GTrackLocus& loctrack, int qcount, GSeqTrack& gtrack, char strand,
    GList<GXLocus>* retxloci=NULL) {
  GList<GXLocus>* dest_xloci=NULL;
  GList<GXLocus> tmpxloci(true,false,true); //local set of newly created xloci
  GList<GXLocus>* xloci=&tmpxloci;
  if (strand=='+') {
       dest_xloci=& gtrack.xloci_f;
       }
    else if (strand=='-') {
      dest_xloci = & gtrack.xloci_r;
      }
   else dest_xloci= & gtrack.xloci_u;

  if (retxloci==NULL) {
     //if no return set of build xloci was given
     //take it as a directive to work directly on the global xloci
     xloci=dest_xloci;
     dest_xloci=NULL;
   }
 for (int q=-1;q<qcount;q++) {
   GList<GLocus>* wrkloci=NULL;
   if (q<0) {
      if (loctrack.rloci.Count()==0) continue;
      //loci=new GList<GLocus>(true,false,false);
      //loci->Add(loctrack.rloc);
      wrkloci = &(loctrack.rloci);
      }
     else {
      if (loctrack[q]==NULL) continue;
      wrkloci = &(loctrack[q]->qloci);
      }

   for (int t=0;t<wrkloci->Count();t++) {
      GLocus* loc=wrkloci->Get(t);
      int xfound=0; //count of parent xloci
      if (loc->xlocus!=NULL) continue; //already assigned a superlocus
      GArray<int> mrgxloci(true);
      for (int xl=0;xl<xloci->Count();xl++) {
         GXLocus& xloc=*(xloci->Get(xl));
         if (xloc.start>loc->end) {
            if (xloc.start-loc->end > GFF_MAX_LOCUS) break;
            continue;
            }
         if (loc->start>xloc.end) continue;
         if (xloc.add_Locus(loc)) {
            xfound++;
            mrgxloci.Add(xl);
            }
         } //for each existing Xlocus
      if (xfound==0) {
         xloci->Add(new GXLocus(loc));
         }
      else {
         int il=mrgxloci[0];
         GXLocus& xloc=*(xloci->Get(il));
         if (xfound>1) {
            for (int l=1;l<xfound;l++) {
              int mlidx=mrgxloci[l]-l+1;
              xloc.addMerge(*(xloci->Get(mlidx)));
              GXLocus* ldel=xloci->Get(mlidx);
              xloci->Delete(mlidx);
              if (retxloci!=NULL)
                    delete ldel;
              }
            }
         //in case xloc.start was decreased, bubble-down until it's in the proper order
         while (il>0 && xloc<*(xloci->Get(il-1))) {
            il--;
            xloci->Swap(il,il+1);
            }
         } //at least one locus is being merged
      }//for each locus
   }//for each set of loci in the region (refs and each qry set)
  //-- add xloci to the global set of xloci unless retxloci was given,
  if (retxloci!=NULL) retxloci->Add(*xloci);
                 else dest_xloci->Add(*xloci);
}

void singleQData(GList<GLocus>& qloci, GList<GTrackLocus>& loctracks) {
 for (int i=0;i<qloci.Count();i++) {
  if (qloci[i]->t_ptr==NULL) {
    GTrackLocus* tloc=new GTrackLocus(numQryFiles);
    tloc->addQLocus(qloci[i],0);
    loctracks.Add(tloc);
    }
  }
}
/*
void recheckUmrnas(GSeqData* gseqdata, GList<GffObj>& mrnas,
     GList<GLocus>& loci, GList<GLocus>& nloci,  GList<GLocus>& oloci) {
 GList<GLocus> reassignedLocs(false,false);
 for (int u=0;u<gseqdata->umrnas.Count();u++) {
   for (int l=0;l<oloci.Count();l++) {
     if (gseqdata->umrnas[u]==NULL) break;
     if (gseqdata->umrnas[u]->end<oloci[l]->start) break; //try next umrna
     if (oloci[l]->end<gseqdata->umrnas[u]->start) continue; //try next locus
     if (gseqdata->umrnas[u]->strand=='+' || gseqdata->umrnas[u]->strand=='-') {
       gseqdata->umrnas.Forget(u);
       continue; //already reassigned earlier
       }
     //umrna overlaps locus region
     GffObj* umrna=gseqdata->umrnas[u];
     for (int m=0;m<oloci[l]->mrnas.Count();m++) {
        if (oloci[l]->mrnas[m]->exonOverlap(umrna)) {
            gseqdata->umrnas.Forget(u);
            CTData* umdata=((CTData*)umrna->uptr);
            //must be in a Loci anyway
            if (umdata==NULL || umdata->locus==NULL)
                GError("Error: no locus pointer for umrna %s!\n",umrna->getID());
            for (int i=0;i<umdata->locus->mrnas.Count();i++) {
               GffObj* um=umdata->locus->mrnas[i];
               um->strand=oloci[l]->mrnas[m]->strand;
               }
            reassignedLocs.Add(umdata->locus);
            break;
            }
        } //for each mrna in locus
      } //for each locus
   } //for each umrna
 if (reassignedLocs.Count()>0) {
   gseqdata->umrnas.Pack();
   gseqdata->nloci_u.setFreeItem(false);
   for (int i=0;i<reassignedLocs.Count();i++) {
     GLocus* loc=reassignedLocs[i];
     for (int m=0;m<loc->mrnas.Count();m++) {
        mrnas.Add(loc->mrnas[m]);
        }
     loci.Add(loc);
     nloci.Add(loc);
     gseqdata->nloci_u.Remove(loc);
     }
   gseqdata->nloci_u.setFreeItem(true);
   }
}
*/

void umrnasXStrand(GList<GXLocus>& xloci, GSeqTrack& gtrack) {
  //try to assign a strand to unoriented transfrags due to overlaps
  //with oriented loci from other samples
 for (int x=0;x<xloci.Count();x++) {
   char newStrand=xloci[x]->strand;
   if (newStrand=='.') continue;
   if (xloci[x]->qloci.Count()==0) continue;
   //go through all qloci in this xlocus
   for (int l = 0; l < xloci[x]->qloci.Count(); l++) {
     char locstrand=xloci[x]->qloci[l]->mrna_maxcov->strand;
     if (locstrand=='.') {
        //this is a unoriented cluster which is now being assigned to xloci[x]->strand
        GLocus* qloc=xloci[x]->qloci[l]; //qry cluster having all its transfrags assigned a strand
        GSeqData* qtrackdata=gtrack.qdata[qloc->qfidx];
        //all the previously unoriented transfrags in qloc
        //will be taken out of qtrackdata->umrnas,
        // but NOTE: their loci are still left in qtrackdata->nloci_u
        // and qtrackdata->loci_f/r are NOT updated
        for (int i=0;i<qloc->mrnas.Count();i++) {
           qloc->mrnas[i]->strand=newStrand; //assign new strand and move
           int uidx=qtrackdata->umrnas.IndexOf(qloc->mrnas[i]);
           if (uidx>=0) {
        	   qtrackdata->umrnas.Forget(uidx);
        	   qtrackdata->umrnas.Delete(uidx);
                if (xloci[x]->strand=='+')
                	qtrackdata->mrnas_f.Add(qloc->mrnas[i]);
                   else
                	qtrackdata->mrnas_r.Add(qloc->mrnas[i]);
                }
           }
        } //unknown strand
     } //for each xloci[x].qloci (l)

   } //for each xloci (x)
}

//cluster loci across all datasets
void xclusterLoci(int qcount, char strand, GSeqTrack& gtrack) {
  //gtrack holds data for all input qry datasets for a chromosome/contig
  //cluster QLoci
  GList<GTrackLocus> loctracks(true,true,false);
  //all vs all clustering across all qry data sets + ref
  //one-strand set of loci from all datasets + ref loci
  GList<GLocus>* wrkloci=NULL;
  //build xloci without references first
  //then add references only if they overlap an existing xloci

  int nq=0;
  for (int q=0;q<=qcount+1;q++) {
    bool refcheck=false;
    if (q==qcount) { // check the unoriented loci for each query file
       while (nq<qcount &&
              (gtrack.qdata[nq]==NULL || gtrack.qdata[nq]->nloci_u.Count()==0))
                 nq++; //skip query files with no unoriented loci
       if (nq<qcount) {
             wrkloci=&(gtrack.qdata[nq]->nloci_u);
             nq++;
             if (nq<qcount) q--; //so we can fetch the next nq in the next q cycle
             }
          else continue; //no more q files with unoriented loci
       }
    else if (q==qcount+1) { // check the reference loci
           if (strand=='+') wrkloci=gtrack.rloci_f;
                       else wrkloci=gtrack.rloci_r;

           if (wrkloci==NULL) break; //no ref loci here
           refcheck=true;
           }
     else  {
          if (gtrack.qdata[q]==NULL) continue;
          if (strand=='+') wrkloci=&(gtrack.qdata[q]->loci_f);
                      else wrkloci=&(gtrack.qdata[q]->loci_r);
         }
   // now do the all-vs-all clustering thing:
   for (int t=0;t<wrkloci->Count();t++) {
      GLocus* loc=wrkloci->Get(t);
      int xfound=0; //count of parent loctracks
      if (loc->t_ptr!=NULL) continue; //already assigned a loctrack
      GArray<int> mrgloctracks(true);
      for (int xl=0;xl<loctracks.Count();xl++) {
         GTrackLocus& trackloc=*loctracks[xl];
         if (trackloc.start>loc->end) break;
         if (loc->start>trackloc.end) continue;
         if (trackloc.add_Locus(loc)) {
            xfound++;
            mrgloctracks.Add(xl);
            }
         } //for each existing Xlocus
      if (xfound==0) {
         if (!refcheck) //we really don't care about ref-only clusters
           loctracks.Add(new GTrackLocus(numQryFiles, loc));
         }
      else {
         int il=mrgloctracks[0];
         GTrackLocus& tloc=*(loctracks.Get(il));
         if (xfound>1) {
           for (int l=1;l<xfound;l++) {
             int mlidx=mrgloctracks[l]-l+1;
             tloc.addMerge(loctracks[mlidx], qcount, loc);
             loctracks.Delete(mlidx);
             }
           }
         //in case tloc.start was decreased, bubble-down 'til it's in the proper place
         while (il>0 && tloc<*(loctracks[il-1])) {
            il--;
            loctracks.Swap(il,il+1);
            }
        } //at least one locus found
      }//for each wrklocus
     } //for each set of loci (q)
   //loctracks is now set with all x-clusters on this strand
 for (int i=0;i<loctracks.Count();i++) {
   if (!loctracks[i]->hasQloci) continue; //we really don't care here about reference-only clusters
   GTrackLocus& loctrack=*loctracks[i];
   findTMatches(loctrack, qcount); //find matching transfrags in this xcluster
   for (int rl=0; rl < loctrack.rloci.Count(); rl++) {
      findTRMatch(loctrack, qcount, *(loctrack.rloci[rl]));
      //find matching reference annotation for this xcluster and assign class codes to transfrags
   }
   GList<GXLocus> xloci(false,false,false);
   buildXLoci(loctrack, qcount, gtrack, strand, &xloci);
   //the newly created xloci are in xloci
   umrnasXStrand(xloci, gtrack);
   //also merge these xloci into the global list of xloci
   for (int l=0; l < xloci.Count(); l++) {
       if (xloci[l]->strand=='+') {
           gtrack.xloci_f.Add(xloci[l]);
           }
          else if (xloci[l]->strand=='-') {
              gtrack.xloci_r.Add(xloci[l]);
              }
            else gtrack.xloci_u.Add(xloci[l]);
   }
 }//for each xcluster
}

//writing .refmap file
void printRefMap(FILE** frs, int qcount, GList<GLocus>* rloci) {
  if (rloci==NULL) return;

  for (int l=0;l<rloci->Count(); l++) {
    for (int r=0;r<rloci->Get(l)->mrnas.Count(); r++) {
      GffObj& ref = *(rloci->Get(l)->mrnas[r]);
      CTData* refdata = ((CTData*)ref.uptr);
      GStr* clist = new GStr[qcount];
      GStr* eqlist = new GStr[qcount];
      for (int i = 0; i<refdata->ovls.Count(); i++) {
        GffObj* m=refdata->ovls[i]->mrna;
        char ovlcode=refdata->ovls[i]->code;
        if (m==NULL) {
          GMessage("Warning: NULL mRNA found for ref %s with ovlcode '%c'\n",
               ref.getID(), refdata->ovls[i]->code);
          continue;
        }
        int qfidx = ((CTData*)m->uptr)->qset;
        if (ovlcode == '=') {
          eqlist[qfidx].append(getGeneID(m));
          eqlist[qfidx].append('|');
          eqlist[qfidx].append(m->getID());
          eqlist[qfidx].append(',');
        }
        else if (ovlcode == 'c') {
          clist[qfidx].append(getGeneID(m));
          clist[qfidx].append('|');
          clist[qfidx].append(m->getID());
          clist[qfidx].append(',');
        }
      }//for each reference overlap
      for (int q=0;q<qcount;q++) {
        if (!eqlist[q].is_empty()) {
          eqlist[q].trimR(',');
          fprintf(frs[q],"%s\t%s\t=\t%s\n", getGeneNameID(ref), ref.getID(),eqlist[q].chars());
        }
        if (!clist[q].is_empty()) {
          clist[q].trimR(',');
          fprintf(frs[q],"%s\t%s\tc\t%s\n",getGeneNameID(ref), ref.getID(),clist[q].chars());
        }
      }
      delete[] clist;
      delete[] eqlist;
    }// ref loop
  }//ref locus loop
}

void trackGData(int qcount, GList<GSeqTrack>& gtracks, GStr& fbasename, FILE** ftr, FILE** frs) {
  FILE* f_ltrack=NULL;
  FILE* f_itrack=NULL;
  FILE* f_ctrack=NULL;
  FILE* f_xloci=NULL;
  int cnum=0; //consensus numbering for printITrack()
  GStr s=fbasename;
  //if (qcount>1 || generic_GFF) { //doesn't make much sense for only 1 query file
    s.append(".tracking");
    f_itrack=fopen(s.chars(),"w");
    if (f_itrack==NULL) GError("Error creating file %s !\n",s.chars());
  //  }
  f_ctrack=fopen(consGTF.chars(),"w");
  if (f_ctrack==NULL) GError("Error creating file %s !\n",s.chars());

  FILE* fredundant=NULL;
  if (discardContained) {
	  s=fbasename;
	  s.append(".redundant.gtf");
	  fredundant=fopen(s.chars(),"w");
  }
  s=fbasename;
  s.append(".loci");
  f_xloci=fopen(s.chars(),"w");
  if (f_xloci==NULL) GError("Error creating file %s !\n",s.chars());
  for (int g=0;g<gtracks.Count();g++) { //for each genomic sequence
    GSeqTrack& gseqtrack=*gtracks[g];

    xclusterLoci(qcount,  '+', gseqtrack);
    xclusterLoci(qcount,  '-', gseqtrack);

    //count XLoci, setting their id
    numXLoci(gseqtrack.xloci_f, xlocnum);
    numXLoci(gseqtrack.xloci_r, xlocnum);
    numXLoci(gseqtrack.xloci_u, xlocnum);
    //transcript accounting: for all those transcripts with 'u' or 0 class code
    // we have to check for polymerase runs 'p' or repeats 'r'

    GFaSeqGet *faseq=gfasta.fetch(gseqtrack.get_gseqid(), checkFasta);

    umrnaReclass(qcount, gseqtrack, ftr, faseq);

    // print transcript tracking (ichain_tracking)
    //if (qcount>1)
    for (int q=0;q<qcount;q++) {
         if (gseqtrack.qdata[q]==NULL) continue;
         printITrack(f_itrack, gseqtrack.qdata[q]->mrnas_f, qcount, cnum);
         printITrack(f_itrack, gseqtrack.qdata[q]->mrnas_r, qcount, cnum);
         //just for the sake of completion:
         printITrack(f_itrack, gseqtrack.qdata[q]->umrnas, qcount, cnum);
         }
    //print XLoci and XConsensi within each xlocus
    //also TSS clustering assignment for XConsensi
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_f, fredundant); //faseq, fredundant);
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_r, fredundant); //faseq, fredundant);
    printXLoci(f_xloci, f_ctrack, qcount, gseqtrack.xloci_u, fredundant); //faseq, fredundant);
    if (tmapFiles && haveRefs) {
      printRefMap(frs, qcount, gseqtrack.rloci_f);
      printRefMap(frs, qcount, gseqtrack.rloci_r);
      }
    delete faseq;
    }
  if (tmapFiles) {
   for (int q=0;q<qcount;q++) {
        fclose(ftr[q]);
        if (haveRefs) fclose(frs[q]);
        }
   }
  if (f_ltrack!=NULL) fclose(f_ltrack);
  if (f_itrack!=NULL) fclose(f_itrack);
  if (f_ctrack!=NULL) fclose(f_ctrack);
  if (f_xloci!=NULL) fclose(f_xloci);

}
