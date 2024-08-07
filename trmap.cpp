#include "GArgs.h"
#include "gff.h"
#include "GStr.h"
#include "GBitVec.h"
#include "GIntervalTree.hh"

#define VERSION "0.12.9"
#define MIN_GFF_VERSION 129

#ifndef GFF_VERSION
 #define GFF_VERSION 0
#endif
#if GFF_VERSION < MIN_GFF_VERSION
 #error "gff.h version mismatch! Please pull/clone gclib"
#endif


bool simpleOvl=false;
bool stricterMatching=false;
bool cdsMatching=false;
bool showCDS=false;
bool outRefOvlTab=false;
bool novelJTab=false;
bool tbest=false;
bool selfMap=false; // --self

GStr fltCodes;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: +strand, 2: -strand
};

const char* USAGE =
"trmap v" VERSION " : transcript to reference mapping and overlap classifier.\nUsage:\n"
"  trmap [-c 'codes'] [-T | -J | -S] <ref_gff> [<query_gff>|--self] \n"
"                    [-B] [-t <Tfile>] [-o <outfile>] "
"\nPositional arguments:\n"
"  <ref_gff>    reference annotation file name (GFF/BED format)\n"
"  <query_gff>  query file name (GFF/BED format) or \"-\" for stdin\n"
"               can be omitted if --self is provided\n"
"Options:\n"
"  -o <outfile> write output to <outfile> instead of stdout\n"
"  --show-cds   add CDS:start:end info to all transcripts with CDS\n"
"  --strict-match : '=' overlap code is assigned when all exons match,\n"
"               while '~' code is assigned when only introns match\n"
"  -c '<codes>' only show overlaps with code in '<codes>' (e.g. -c '=ck')\n"
"  -T           output a 7 column table with overlap info: \n"
"                 queryID, class_code, refID, ref_cov, rev_ovl_bias,\n"
"                    matching_ref_introns, num_matching_junctions\n"
"                ..where matching_ref_introns has this format:\n"
"                   total_ref_intron_count:list_of_matching_introns\n"
"  --best, -B   for -T/-t option, only show the \"best\" class code\n"
"  -t <Tfile>   write the table described for -T to file <Tfile>, while\n"
"               allowing other option for main output\n"
"  -J           for each query transcript output a 6 column table\n"
"                 queryID, chr:strand:exons, list of reference transcripts,\n"
"                 num ref genes, list of gene names, list of novel junctions\n"
"  -S           report only simple exon overlap percentage with any reference\n"
"               transcripts (one line per query)\n";

//"  -f <file>    for -J, report \"fusion\" query transcripts, i.e. transcripts\n"
//"               that overlaps multiple non-overlapping reference genes \n"

bool closerRef(GffObj* a, GffObj* b, int numexons, byte rank) {
 //this is called when a query overlaps a and b with the same overlap length
 //to decide which of a or b is closer structurally to the query
 // returns true if a is closer, false if b is closer
 if (a==NULL || b==NULL) return (a!=NULL);
 if (rank<CLASSCODE_OVL_RANK) {
	 //significant intron/exon overlap -- all the 'j' codes, but includes 'e'
	 if (a->exons.Count()!=b->exons.Count()) {
		 int ad=a->exons.Count()-numexons;
		 int bd=b->exons.Count()-numexons;
		 return (abs(ad)==abs(bd)) ? ad<bd : abs(ad) < abs(bd);
	 }
 }
 if (a->exons.Count()!=b->exons.Count()) return (a->exons.Count()>b->exons.Count());
 if (a->hasCDS() && !b->hasCDS())
        return true;
   else {
     if (b->hasCDS() && !a->hasCDS()) return false;
     return (a->covlen==b->covlen) ? (strcmp(a->getID(), b->getID())<0) :
    		 (a->covlen>b->covlen);
     }
 }


struct TRefOvl {
	GffObj* ref;
    byte rank;
	//char ovlcode;
    //int ovlen;
    TOvlData* od; //copy of overlap data
    int numExons; //number of exons in the query mRNA

    bool operator<(TRefOvl& b) { //lower = higher priority
		if (rank==b.rank && b.ref!=NULL && ref!=NULL) {
			if (numExons==1 && b.ref->exons.Count()!=ref->exons.Count()) {
				if (ref->exons.Count()==1) return true; //SET match always has priority
				else if (b.ref->exons.Count()==1) return false;
			}
			if (numExons>1 && b.od->numJmatch!=od->numJmatch) {
				return (od->numJmatch > b.od->numJmatch);
			}
			return (od->ovlen==b.od->ovlen) ?
				closerRef(ref, b.ref, numExons, rank) : (od->ovlen>b.od->ovlen);
		}
		else return rank<b.rank;
	}
    bool operator==(TRefOvl& b) {
		return (rank==b.rank && ref==b.ref);
	}

	TRefOvl(GffObj* r, TOvlData& o, int exonCount=0): ref(r),
		rank(classcode_rank(o.ovlcode)), numExons(exonCount) {
		od=new TOvlData(o);
	}
	~TRefOvl() { delete od; }
};

struct QJData {
	GBitVec jmd; //bit array for junctions (1 = ref-matched junction)
	GBitVec inmd; //bit array for introns (1 = ref-matched intron)
    //exon-skipping events = 11 in jmd having no corresponding 1 bit set in inmd
	GList<TRefOvl> refovls;
	GffObj* t; //query transcript
	GVec<GSeg> introns;
	QJData():t(NULL) {}
	QJData(GffObj& tr):jmd( (tr.exons.Count()-1)<<1 ), inmd(tr.exons.Count()-1),
			refovls(true,true,true), t(&tr), introns() {
		GSeg in(0,0);
		for (int i=1;i<tr.exons.Count();i++) {
			in.start=tr.exons[i-1]->end+1;
			in.end=tr.exons[i]->start-1;
			introns.Add(&in);
		}
	}
	// return index of intron matching given intron coordinates if any
	// assumes that introns are sorted, non-overlapping
    int findIntron(uint istart, uint iend, int i0=0) {
    	int r=-1;
    	for (int i=i0;i<introns.Count();i++) {
    		if (introns[i].start>istart) break;
    		if (introns[i].start==istart && introns[i].end==iend)
    		   { r=i; break; }
    	}
    	return r;
    }

	void add(GffObj* ref, TOvlData& od) {
	    #ifndef NDEBUG
		  if (jmd.size()!=od.jbits.size()) GError("Error: mismatching QJData bit vector size!\n");
		#endif
		bool sameStrand=(t->strand == ref->strand);
		if (sameStrand) {
		   jmd |= od.jbits;
           inmd |= od.inbits;
		}
		int idx=refovls.Add(new TRefOvl(ref, od, t->exons.Count()));
		#ifndef NDEBUG
		  if (idx<0) {
			  refovls.Found(new TRefOvl(ref, od, t->exons.Count()), idx);
			  GMessage("Error: %s trying to add duplicate overlap of ref %s (previously found as %s)!\n",
					  t->getID(), ref->getID(), refovls[idx]->ref->getID() );
			  GMessage(" Existing ref overlaps:\n");
			  for (int i=0;i<refovls.Count();i++) {
				  GMessage("%c\t%s\n", refovls[i]->od->ovlcode, refovls[i]->ref->getID());
			  }
			  GError("exiting..\n");
		  }
		#endif
	}
};

struct RefGene:public GSeg { //start/end: min-max observed for transcripts in this gene
  const char* gene_id; //points to gene_ids in ref GffObj
  // keep track of ref genes with their merged exons and ranges
  GArray<GSeg> mexons; //set of non-overlapping exon spans (resulted after merging overlapping exons)
  GVec<GffObj*> transcripts;
  //GVec<GSeg> introns; //set of all introns (could be overlapping)
  GVec<RefGene*> ogenes; //overlapping genes (based on mexon-overlaps only)
  RefGene():gene_id(NULL),mexons(true, true),transcripts(), ogenes() {}
  RefGene(const char* gid, GffObj& t):GSeg(t.start, t.end),
		  gene_id(gid),mexons(true, true), transcripts(), ogenes() {
	  transcripts.cAdd(&t);
	  for (int i=0;i<t.exons.Count();i++)  {
			  mexons.Add(t.exons[i]);
	  }
  }

  int getLength() {
	int r=0;
	for (int i=0;i<mexons.Count();i++)
		r+=mexons[i].len();
	return r;
  }

  void addTranscript(GffObj& t) {
	  transcripts.cAdd(&t);
	  for (int i=0;i<t.exons.Count();i++) { //for each exon of the incoming transcript
		  int ni=mexons.Add(t.exons[i]); //exon is added to mexons as ordered by start coord
		  if (ni>=0) { //exon added, check overlaps with other mexons and merge as needed
			  int delcount=0;
			  int j=ni+1;
			  //could overlap the preceding mexon (mexons[ni-1])
			  if (ni>0 && mexons[ni-1].end>=mexons[ni].start) {
				  ni--; //merge into preceding exon, which now takes over
				  delcount++;
			  }
			  //check how many mexons are being overlapped by this new added exon
			  while (j<mexons.Count() && mexons[ni].end>mexons[j].start) {
				 delcount++;
                 j++;
			  }
			  if (delcount) {
				  j--; //merge with last overlapped range
				  if (mexons[ni].end<mexons[j].end)
					   mexons[ni].end=mexons[j].end;
				  mexons.Delete(ni+1,delcount);
			  }
		  }
	  }
  }
};

//queryID, chr:strand:exons, list of ovlcode|ref_transcripts|gene,
//                num genes,  list of novel junctions

void geneAdd(GVec<char*>& gset, char* g) {
	for (int i=0;i<gset.Count();i++) {
		if (strcmp(g, gset[i])==0) return;
	}
	gset.Add(g);

}

// for sorting GVec(<char*>):
int cmpcstr(const pointer p1, const pointer p2) {
 return strcmp(*(char**)p1, *(char**)p2);
}

void printNJTab(FILE* f, QJData& d) {
	fprintf(f, "%s\t%s:%c", d.t->getID(), d.t->getRefName(), d.t->strand);
	for (int i=0;i<d.t->exons.Count();++i) {
		char ch=i ? ',' : ':';
		fprintf(f, "%c%d-%d", ch, d.t->exons[i]->start, d.t->exons[i]->end);
	}

	fprintf(f, "\t");

	if (d.refovls.Count()==0) { //no ref range overlaps found
		const char* code=d.t->getAttr("class_code"); //preserve u,p,r info
		if (code==NULL) code=".";
		fprintf(f,"%s\t.\t",code);
		if (d.t->exons.Count()<=1) fprintf(f, ".");
		else for (int i=1;i<d.t->exons.Count();i++) {
			   //every junction is going to be novel:
			   if (i>1) fprintf(f, ",");
			   fprintf(f, "%d-%d:nn", d.t->exons[i-1]->end+1, d.t->exons[i]->start-1);
			 }
		fprintf(f, "\n");
		return;
	}
	GVec<char*> genes; //gene IDs
	// for self mapping, add own gene name
	if (selfMap) geneAdd(genes, d.t->getGeneName());
	for (int i=0;i<d.refovls.Count();++i) {
		char* g=d.refovls[i]->ref->getGeneName();
		if (g==NULL) g=d.refovls[i]->ref->getGeneID();
		if (i) fprintf(f, ",");
		fprintf(f, "%c|%s|", d.refovls[i]->od->ovlcode, d.refovls[i]->ref->getID());
		if (g) {
			fprintf(f, "%s", g);
			if (d.refovls[i]->rank<CLASSCODE_OVL_RANK)
			    geneAdd(genes, g); //only count actually overlapping
		} else fprintf(f, ".");
	}
	//fprintf(f, "\t%d\t", genes.Count());
	fprintf(f, "\t");
	//print list of gene names if possible
	if (genes.Count()==0) fprintf(f, ".");
	else {
		//alpha sort if multiple
		if (genes.Count()>1)
			genes.Sort(&cmpcstr);
		for (int i=0;i<genes.Count();i++) {
			if (i) fprintf(f, ",");
			fprintf(f, "%s", genes[i]);
		}
	}
	fprintf(f, "\t");
	// now print novel junctions, in groups of 2:nn|n.|.n|ss|..
	//ss is exon skip code = novel intron, even though both splice sites are known
	char jj[3]={'.','.','\0'};
	bool printed=false;
	for (uint i=0;i<d.jmd.size();i+=2) {
		bool smatch=d.jmd[i];
		bool ematch=d.jmd[i+1];
		if (smatch && ematch) {
			if (d.inmd[i>>1]) continue; //known intron
			//exon skipping! novel intron without novel junctions
			jj[0]= 's';
			jj[1]= 's';
		} else {
			jj[0]= (smatch) ? '.' : 'n';
			jj[1]= (ematch) ? '.' : 'n';
		}
		int ei = i>>1; // index of exon on the left
		if (printed) fprintf(f, ",");
		printed=true;
		fprintf(f, "%d-%d:%s", d.t->exons[ei]->end+1, d.t->exons[ei+1]->start-1, jj);
	}
	if (!printed) fprintf(f, ".");
	fprintf(f, "\n");

}

void printOvlTab(FILE* fwtab, const char* tid, GffObj* r, TOvlData& od, const char* tgene) {
	if (!r) return;
	const char* rgi=r->getGeneID(); if (rgi==NULL) rgi="";
	const char* rgn=r->getGeneName(); if (rgn==NULL) rgn="";
	if (tgene==NULL) tgene="";
	fprintf(fwtab, "%s|%s\t%c\t%s|%s|%s", tid, tgene, od.ovlcode, r->getID(), rgi ,rgn );
	if (od.ovlen) {
		float rcov= (100.00*od.ovlen)/r->covlen;
		fprintf(fwtab, "\t%1.f\t", rcov);
	} else fprintf(fwtab, "\t.\t");
    float rovlbias=0;
	if (r->strand=='-') {
		if (od.ovlen) {
			int rs=r->covlen-od.ovlen-(od.ovlRefstart-1);
			int rsmid=rs+od.ovlen/2;
			int rmid=r->covlen/2;
			//overlap deviation from center (centered=0.5)
			rovlbias=0.5+((float)(rsmid-rmid))/r->covlen;
			fprintf(fwtab, "%.2f", rovlbias);
		} else fprintf(fwtab,".");
		int im=-1;
		int nint=od.rint.size(); //number of introns
		if (nint)
			im=od.rint.find_first(); //find first matching intron in ref
		if (im>=0) {
		   GVec<int> introns;
		   introns.cAdd(nint-im);
		   while( (im=od.rint.find_next(im))>0 )
			   introns.cAdd(nint-im);
		   im=introns.Pop();
		   fprintf(fwtab,"\t%d:%d", nint, im );
		   int ni=introns.Count()-1;
		   for (int i=ni;i>=0;--i)
			   fprintf(fwtab,",%d", introns[i]);
		} else fprintf(fwtab, "\t.");
	} else {
		if (od.ovlen) {
			int rsmid=od.ovlRefstart-1+od.ovlen/2;
			int rmid=r->covlen/2;
			//overlap deviation from center (centered=0.5)
			rovlbias=0.5+((float)(rsmid-rmid))/r->covlen;
			fprintf(fwtab, "%.2f", rovlbias);
		} else fprintf(fwtab,".");
		int im=-1;
		int nint=od.rint.size();
		if (nint)
			im=od.rint.find_first();
		if (im>=0) {
			   fprintf(fwtab,"\t%d:%d", nint, im+1 );
		   while( (im=od.rint.find_next(im))>0 )
			   fprintf(fwtab,",%d", im+1);
		} else fprintf(fwtab, "\t.");
	}
	fprintf(fwtab, "\t%d\n", od.numJmatch);
}

void printTabBest(FILE* fwtab, QJData& d) {
  if (d.refovls.Count()>0) {
	  TRefOvl* ro=d.refovls.First();
	  printOvlTab(fwtab, d.t->getID(), ro->ref, *(ro->od), d.t->getGeneName());
  }
}

int nextQi=0;
GffObj* getNextQtx(GffReader* myQ, GPVec<GffObj> *refKeep) {
  GffObj* t=NULL;
  if (myQ) return myQ->readNext();
  if (nextQi<refKeep->Count()) {
	  t= refKeep->Get(nextQi);
	  ++nextQi;
  }
  return t;
}

int main(int argc, char* argv[]) {
	FILE* fwtab=NULL;
	GArgs args(argc, argv, "help;strict-match;best;show-cds;self;hTBJSc:t:o:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help")) {
		GMessage(USAGE);
		exit(EXIT_SUCCESS);
	}
	bool optT=false; //-T output ONLY
	if (args.getOpt('S')) simpleOvl=true;
	if (args.getOpt("strict-match")) stricterMatching=true;
	if (args.getOpt("show-cds")) showCDS=true;
	if (args.getOpt("self")) selfMap=true;
	if (args.getOpt("best") || args.getOpt('B')) tbest=true;
	if (args.getOpt('T')) outRefOvlTab=optT=true;
	if (args.getOpt('J')) novelJTab=true;
	if ((int)outRefOvlTab + (int)simpleOvl+(int)novelJTab > 1)
		GError("%s\nError: options -T, -J and -S are mutually exclusive!\n", USAGE);
	const char* tfo=args.getOpt('t');
	if (tfo) { //secondary -t output
		outRefOvlTab=true;
		if (strcmp(tfo, "-")==0) fwtab=stdout;
			          else {
			            fwtab=fopen(tfo, "w");
			            if (fwtab==NULL) GError("Error creating output file %s !\n",tfo);
			          }
	}
	if (tbest && !outRefOvlTab)
		GError("%s\nError: option -B/--best requires -T or -t!\n", USAGE);
	// tbest implies outRefOvlTab	
    const char*s=args.getOpt('c');
    if (s!=NULL) {
    	fltCodes=s;
    }
	GHash<GSTree*> map_trees; //map a ref sequence name to its own interval trees (3 per ref seq)

	const char* o_file = args.getOpt('o') ? args.getOpt('o') : "-";
    int pospar=args.startNonOpt();
    if (selfMap) {
    	if (pospar>1) GError("%s\nError: only one transcript file expected with --self!\n");
    	if (pospar<1) GError("%s\nError: a transcript file is expected as input!\n");
    } else if (pospar!=2)
		GError("%s\nError: %d arguments provided (expected 2)\n",USAGE, args.startNonOpt());
	const char* ref_file = args.nextNonOpt();
	const char* q_file = args.nextNonOpt();

	FILE* fr=fopen(ref_file, "r");
	if (fr==NULL) GError("Error: could not open reference annotation file (%s)!\n", ref_file);

	GffReader* myR=new GffReader(fr, true, true);
	if (selfMap && (novelJTab || tbest)) myR->keepAttrs();
	const char* fext=getFileExt(ref_file);
	if (Gstricmp(fext, "bed")==0) myR->isBED();
	GffObj* t=NULL;
	GPVec<GffObj> *refKeep = new GPVec<GffObj>(true);
	while ((t=myR->readNext())!=NULL) {
		if (t->exons.Count()==0) {
			delete t;
			continue; //skip exonless entities (e.g. genes)
		}
		GSTree* cTree=map_trees[t->getGSeqName()];
		if (cTree==NULL) {
			cTree=new GSTree();
			map_trees.Add(t->getGSeqName(), cTree);
		}
		if (t->strand=='+')
		 cTree->it[1].Insert(t);
		else if (t->strand=='-')
			cTree->it[2].Insert(t);
		else cTree->it[0].Insert(t);
		refKeep->Add(t);
	}
	delete myR;
	FILE* outFH=NULL;
	if (strcmp(o_file, "-")==0) outFH=stdout;
	            else {
	            	outFH=fopen(o_file, "w");
	            	if (outFH==NULL) GError("Error creating file %s !\n",o_file);
	            }
	if (outRefOvlTab && fwtab==NULL) fwtab=outFH;
	FILE* fq=NULL;
	fext=NULL;
	GffReader* myQ = NULL;
	if (q_file!=NULL) {
		if (strcmp(q_file,"-")==0) fq=stdin;
		else {
			fq=fopen(q_file, "r");
			if (fq==NULL)
				GError("Error: could not open query file (%s)!\n", q_file);
			fext=getFileExt(q_file);

			myQ = new GffReader(fq, true, true);
			if (novelJTab || tbest) myQ->keepAttrs();
			if (fext && Gstricmp(fext, "bed")==0) myQ->isBED();
		}
	} else if (!selfMap) GError("Error: --self required with only one input file!\n");

	t=NULL;
	while ((t=getNextQtx(myQ, refKeep))!=NULL) {
		const char* gseq=t->getGSeqName();
		if (!map_trees.hasKey(gseq)) {
			if (!selfMap) delete t;
			continue; //reference sequence not present in annotation, so we can't compare
		}
		if (t->exons.Count()==0) {
			if (!selfMap) delete t;
			continue; //only work with properly defined transcripts
		}
		GSTree* cTree=map_trees[gseq];
		//GVec<int> sidx;
		// always search all strands
		QJData* tjd=NULL;
		if (novelJTab || tbest) tjd=new QJData(*t);
		for (int k=0;k<3;++k) {
			GVec<GSeg*> *enu = cTree->it[k].Enumerate(t->start, t->end);
			if(enu->Count()==0) { delete enu; continue; } // no overlaps found
			bool qprinted=false;
			for (int i=0; i<enu->Count(); ++i) { //for each range overlap
				GffObj* r=(GffObj*)enu->Get(i);
				if (selfMap && strcmp(r->getID(), t->getID())==0)
					continue; // skip self matches
				TOvlData od=getOvlData(*t, *r, stricterMatching, 1, cdsMatching);
				// opposite strand non-overlaps should be ignored ?
				bool Xstrand=(t->strand!=r->strand && t->strand!='.' && r->strand!='.');
				if (Xstrand && classcode_rank(od.ovlcode)<CLASSCODE_OVL_RANK) {
					//for trmap, mark these as non-overlaps, with code 'x'
					od.ovlcode='x';
					od.ovlen=0;
				}
				//no real code found (?)
				if (!fltCodes.is_empty() && !fltCodes.contains(od.ovlcode))
					continue;
				//bool novlXstrand = (Xstrand && classcode_rank(od.ovlcode)<classcode_rank('i'));
				//if (novlXstrand) continue;
				// -- two output modes: aggregating (best/sorted), or as-you-go, for each overlap
				// novelJTab and tbest are aggregating
				if (novelJTab || tbest) {
					tjd->add(r, od);
				}

				if (outRefOvlTab && !tbest) {// -T or -t, but not -B/--best
						//if (tbest) tjd->add(r, od);else 
						printOvlTab(fwtab, t->getID(), r, od, t->getGeneName());
				}
				// could be default pseudo-fasta, or simpleOvl
				if (simpleOvl) { //-S output
						if (Xstrand || od.ovlen==0) continue;
						float rcov=(100.00*od.ovlen)/r->covlen;
						if (!qprinted) {
							fprintf(outFH, "%s\t%s:%d-%d|%c", t->getID(), gseq, t->start, t->end, t->strand);
							qprinted=true;
						}
						//append each overlapping referenced to the same line
						fprintf(outFH, "\t%s:%.1f", r->getID(), rcov);
				} // -S output
				else if (!(optT || novelJTab)) { //no -T, -J or -S => default detailed pseudo-FASTA output
						if (!qprinted) {
							fprintf(outFH, ">%s %s:%d-%d %c ", t->getID(), t->getGSeqName(), t->start, t->end, t->strand);
							t->printExonList(outFH);
							if (showCDS && t->hasCDS()) {
							  fprintf(outFH, " CDS:");
							  t->printCDSList(outFH);
							}
							fprintf(outFH, "\n");
							qprinted=true;
						}
						fprintf(outFH, "%c\t", od.ovlcode);
						fprintf(outFH, "%s\t%c\t%d\t%d\t%s\t", r->getGSeqName(), r->strand,
							r->start, r->end, r->getID());
						r->printExonList(outFH);
						if (showCDS && r->hasCDS()) {
						  fprintf(outFH, "\tCDS:");
						  r->printCDSList(outFH);
						}
						fprintf(outFH, "\n");
					} //default pseudo-FASTA output
				//} //no -J
			} //for each range overlap
			if (simpleOvl && qprinted)
			      fprintf(outFH, "\n"); //for simpleOvl all overlaps are on a single line
			delete enu;
		} //for each searchable strand
		if (tjd) {
			if (novelJTab) {
				printNJTab(outFH, *tjd);
			}
			if (tbest) {
				printTabBest(fwtab, *tjd);
			}
			delete tjd;
		}
		if (!selfMap) delete t;
	}
	delete myQ;
    delete refKeep;
	if (outFH!=stdout) fclose(outFH);
	if (fwtab && fwtab!=outFH && fwtab!=stdout) fclose(fwtab);
	return 0;
}
