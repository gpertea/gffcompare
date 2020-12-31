#include "GArgs.h"
#include "gff.h"
#include "GStr.h"
#include "GIntervalTree.h"

#define VERSION "0.12.2"

bool simpleOvl=false;
bool stricterMatching=false;
bool showCDS=false;
bool outTab=false;
GStr fltCodes;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: +strand, 2: -strand
};

const char* USAGE =
"trmap v" VERSION " : transcript to reference mapping and overlap classifier.\nUsage:\n"
"  trmap [-c 'codes'] [-T | -S] [-o <outfile>] <ref_gff> <query_gff>\n"
"Positional arguments:\n"
"  <ref_gff>    reference annotation file name (GFF/BED format)\n"
"  <query_gff>  query file name (GFF/BED format) or \"-\" for stdin\n"
"Options:\n"
"  -o <outfile> write output to <outfile> instead of stdout\n"
"  --show-cds   add CDS:start:end info to all transcripts with CDS\n"
"  --strict-match : '=' overlap code is assigned when all exons match,\n"
"               while '~' code is assigned when only introns match\n"
"  -c '<codes>' only show overlaps with code in '<codes>' (e.g. -c '=ck')\n"
"  -T           output 4 column table: queryID, ovl_code, ref_cov%, refID\n"
"  -S           report only simple exon overlap percentages with reference\n"
"               transcripts, without classification (one line per query)\n";

int main(int argc, char* argv[]) {
	GArgs args(argc, argv, "help;strict-match;show-cds;hTSc:o:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help")) {
		GMessage(USAGE);
		exit(EXIT_SUCCESS);
	}
	if (args.getOpt('S')) simpleOvl=true;
	if (args.getOpt("strict-match")) stricterMatching=true;
	if (args.getOpt("show-cds")) showCDS=true;
	if (args.getOpt("T")) outTab=true;
	if (outTab && simpleOvl)
		GError("%s\nError: options -S and -T are mutually exclusive!\n", USAGE);
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
						float rcov=(100.00*od.ovlen)/r->len();
						fprintf(outFH, "%s\t%c\t%.1f\t%s\n", t->getID(), od.ovlcode, rcov, r->getID());
					} else { //full pseudo-FASTA output
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
