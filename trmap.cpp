#include "GArgs.h"
#include "gff.h"
#include "GStr.h"
#include "GBitVec.h"
#include "GIntervalTree.h"

#define VERSION "0.12.5"

bool simpleOvl=false;
bool stricterMatching=false;
bool showCDS=false;
bool outTab=false;
bool novelJTab=false;
GStr fltCodes;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: +strand, 2: -strand
};

const char* USAGE =
"trmap v" VERSION " : transcript to reference mapping and overlap classifier.\nUsage:\n"
"  trmap [-c 'codes'] [-T | -J | -S] [-o <outfile>] <ref_gff> <query_gff>\n"
"Positional arguments:\n"
"  <ref_gff>    reference annotation file name (GFF/BED format)\n"
"  <query_gff>  query file name (GFF/BED format) or \"-\" for stdin\n"
"Options:\n"
"  -o <outfile> write output to <outfile> instead of stdout\n"
"  --show-cds   add CDS:start:end info to all transcripts with CDS\n"
"  --strict-match : '=' overlap code is assigned when all exons match,\n"
"               while '~' code is assigned when only introns match\n"
"  -c '<codes>' only show overlaps with code in '<codes>' (e.g. -c '=ck')\n"
"  -T           output a 5 column table: \n"
"                 queryID, ovl_code, ref_cov%, refID, matched junctions\n"
"  -J           for each query transcripts, output an 5 column table\n"
"                queryID, chr:strand:exons, list of reference transcripts,\n"
"                num ref genes, list of novel junctions\n"
"  -S           report only simple exon overlap percentages with reference\n"
"               transcripts, without classification (one line per query)\n";

bool closerRef(GffObj* a, GffObj* b, int numexons, byte rank) {
 //this is called when a query overlaps a and b with the same overlap length
 //to decide which of a or b is closer structurally to the query
 // returns true if a is closer, false if b is closer
 if (a==NULL || b==NULL) return (a!=NULL);
 if (rank<CLASSCODE_OVL_RANK) {
	 //significant intron/exon overlap -- all the 'j' codes, but includes 'e'
	 if (a->exons.Count()!=b->exons.Count())
		 return (abs(a->exons.Count()-numexons)<abs(b->exons.Count()-numexons));
 }
 if (a->exons.Count()!=b->exons.Count()) return (a->exons.Count()>b->exons.Count());
 if (a->hasCDS() && !b->hasCDS())
        return true;
   else {
     if (b->hasCDS() && !a->hasCDS()) return false;
     return (a->covlen>b->covlen);
     }
 }


struct TRefOvl {
	GffObj* ref;
	char ovlcode;
    byte rank;
    int ovlen;
    int16_t numExons; //number of exons in the query mRNA
    int16_t numJmatch; //number of matching junctions in this overlap

    bool operator<(TRefOvl& b) { //lower = higher priority
		if (rank==b.rank && b.ref!=NULL && ref!=NULL) {
			if (numExons==1 && b.ref->exons.Count()!=ref->exons.Count()) {
				if (ref->exons.Count()==1) return true; //SET match always has priority
				else if (b.ref->exons.Count()==1) return false;
			}
			if (numExons>1 && b.numJmatch!=numJmatch) {
				return (numJmatch > b.numJmatch);
			}
			return (ovlen==b.ovlen)? closerRef(ref, b.ref, numExons, rank) : (ovlen>b.ovlen);
		}
		else return rank<b.rank;
	}
    bool operator==(TRefOvl& b) {
		return (rank==b.rank && ref==b.ref);
	}

	TRefOvl(GffObj* r=NULL, char code=0, int exonCount=0, int olen=0, int jmatch=0):ref(r),
			ovlcode(code), ovlen(olen), numJmatch(jmatch) {
    	if (exonCount>255) exonCount=255;
    	numExons=exonCount;
		rank=classcode_rank(code);
	}
};

struct QJData {
	GBitVec jmd; //bit array showing ref-matched junctions
	GList<TRefOvl> refovls;
	GffObj* t;
	QJData():t(NULL) {}
	QJData(GffObj& tr):jmd( (tr.exons.Count()-1)<<1 ), refovls(true,true,true),
			t(&tr) { }
	void add(GffObj* ref, TOvlData& od) {
	    #ifndef NDEBUG
		  if (jmd.size()!=od.jbits.size()) GError("Error: mismatching QJData bit vector size!\n");
		#endif
		jmd |= od.jbits;
		int idx=refovls.Add(new TRefOvl(ref, od.ovlcode, t->exons.Count(), od.ovlen, od.numJmatch));
		#ifndef NDEBUG
		  if (idx<0) GError("Error: trying to add duplicate overlap of %s with !\n");
		#endif
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
void printNJTab(FILE* f, QJData& d) {
	if (d.refovls.Count()==0) return;
	fprintf(f, "%s\t%s:%c", d.t->getID(), d.t->getRefName(), d.t->strand);
	GVec<char*> genes; //gene IDs
	for (int i=0;i<d.t->exons.Count();++i) {
		char ch=i ? ',' : ':';
		fprintf(f, "%c%d-%d", ch, d.t->exons[i]->start, d.t->exons[i]->end);
	}
	fprintf(f, "\t");
	for (int i=0;i<d.refovls.Count();++i) {
		char* g=d.refovls[i]->ref->getGeneID();
		if (i) fprintf(f, ",");
		fprintf(f, "%c|%s|", d.refovls[i]->ovlcode, d.refovls[i]->ref->getID());
		if (g) {
			geneAdd(genes, g);
			fprintf(f, "%s", g);
		} else fprintf(f, "_");
	}
	fprintf(f, "\t%d\t", genes.Count());
	// now print novel junctions, in groups of 2 :nn|n.|.n
	char jj[3]={'.','.','\0'};
	bool printed=false;
	for (uint i=0;i<d.jmd.size();i+=2) {
		bool smatch=d.jmd[i];
		bool ematch=d.jmd[i+1];
		if (smatch && ematch) continue;
		int ei = i>>1; // index of exon on the left
		jj[0]= (smatch) ? '.' : 'n';
		jj[1]= (ematch) ? '.' : 'n';
		if (printed) fprintf(f, ",");
		printed=true;
		fprintf(f, "%d-%d:%s", d.t->exons[ei]->end+1, d.t->exons[ei+1]->start-1, jj);
	}
	if (!printed) fprintf(f, ".");
	fprintf(f, "\n");

}

int main(int argc, char* argv[]) {
	GArgs args(argc, argv, "help;strict-match;show-cds;hTJSc:o:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help")) {
		GMessage(USAGE);
		exit(EXIT_SUCCESS);
	}
	if (args.getOpt('S')) simpleOvl=true;
	if (args.getOpt("strict-match")) stricterMatching=true;
	if (args.getOpt("show-cds")) showCDS=true;
	if (args.getOpt('T')) outTab=true;
	if (args.getOpt('J')) novelJTab=true;
	if ((int)outTab + (int)simpleOvl+(int)novelJTab > 1)
		GError("%s\nError: options -T, -J and -S are mutually exclusive!\n", USAGE);
    const char* s=args.getOpt('c');
    if (s!=NULL) {
    	fltCodes=s;
    }
	GHash<GSTree*> map_trees; //map a ref sequence name to its own interval trees (3 per ref seq)

	const char* o_file = args.getOpt('o') ? args.getOpt('o') : "-";

	if (args.startNonOpt()!=2)
		GError("%s\nError: %d arguments provided (expected 2)\n",USAGE, args.startNonOpt());
	const char* ref_file = args.nextNonOpt();
	const char* q_file = args.nextNonOpt();

	FILE* fr=fopen(ref_file, "r");
	if (fr==NULL) GError("Error: could not open reference annotation file (%s)!\n", ref_file);

	GffReader* myR=new GffReader(fr, true, true);
	const char* fext=getFileExt(ref_file);
	if (Gstricmp(fext, "bed")==0) myR->isBED();
	GffObj* t=NULL;
	GPVec<GffObj> *toFree = new GPVec<GffObj>(true);
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
		toFree->Add(t);
	}
	delete myR;
	FILE* outFH=NULL;
	if (strcmp(o_file, "-")==0) outFH=stdout;
	            else {
	            	outFH=fopen(o_file, "w");
	            	if (outFH==NULL) GError("Error creating file %s !\n",o_file);
	            }
	FILE* fq=NULL;
	fext=NULL;
	if (strcmp(q_file,"-")==0) fq=stdin;
	else {
		fq=fopen(q_file, "r");
		if (fq==NULL)
			GError("Error: could not open query file (%s)!\n", q_file);
		fext=getFileExt(q_file);
	}
	GffReader* myQ = new GffReader(fq, true, true);
	if (fext && Gstricmp(fext, "bed")==0) myQ->isBED();
	t=NULL;
	while ((t=myQ->readNext())!=NULL) {
		const char* gseq=t->getGSeqName();
		if (!map_trees.hasKey(gseq)) {
			delete t;
			continue; //reference sequence not present in annotation, so we can't compare
		}
		if (t->exons.Count()==0) {
			delete t;
			continue; //only work with properly defined transcripts
		}
		GVec<int> sidx;
		sidx.cAdd(0); //always search the '.' strand
		if (t->strand=='+') sidx.cAdd(1);
		else if (t->strand=='-') sidx.cAdd(2);
		else { sidx.cAdd(1); sidx.cAdd(2); }
		for (int k=0;k<sidx.Count();++k) {
			GVec<GSeg*> *enu = map_trees[gseq]->it[sidx[k]].Enumerate(t->start, t->end);

			if(enu->Count()>0) { //overlaps found
				QJData* tjd=NULL;
				if (novelJTab) tjd=new QJData(*t);
				bool qprinted=false;
				for (int i=0; i<enu->Count(); ++i) {
					GffObj* r=(GffObj*)enu->Get(i);
					TOvlData od=getOvlData(*t, *r, stricterMatching);
					if (!fltCodes.is_empty() && !fltCodes.contains(od.ovlcode))
						continue;
					if (simpleOvl) {
						if (od.ovlen==0) continue;
						float rcov=(100.00*od.ovlen)/r->covlen;
						if (!qprinted) {
							fprintf(outFH, "%s\t%s:%d-%d|%c", t->getID(), gseq, t->start, t->end, t->strand);
							qprinted=true;
						}
						//append each overlapping referenced to the same line
						fprintf(outFH, "\t%s:%.1f", r->getID(), rcov);
					} else if (outTab) { //3 column output
						//int xovlen=r->exonOverlapLen(*t);
						//GMessage("DEBUG:: Exon overlap: %d (reported by getOvlData: %d)\n", xovlen, od.ovlen);
						float rcov=(100.00*od.ovlen)/r->covlen;
						fprintf(outFH, "%s\t%c\t%.1f\t%s\t%d\n", t->getID(), od.ovlcode,
								rcov, r->getID(), od.numJmatch);
					} else if (novelJTab) {
						tjd->add(r, od);
					}
					else  { //full pseudo-FASTA output
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
					}
				} //for each range overlap
				if (simpleOvl && qprinted)
					fprintf(outFH, "\n"); //for simpleOvl all overlaps are on a single line
				if (novelJTab) {
					printNJTab(outFH, *tjd);
					delete tjd;
				}
			} //has overlaps
			delete enu;
		}
		delete t;
	}
	delete myQ;
    delete toFree;
	fclose(outFH);
	return 0;
}
