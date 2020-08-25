#include "GArgs.h"
#include "gff.h"
#include "GIntervalTree.h"

#define VERSION "0.11.6"

bool simpleOvl=false;
bool stricterMatching=false;
bool showCDS=false;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: +strand, 2: -strand
};

const char* USAGE = "trmap v" VERSION " : transcript to reference mapping and overlap classifier.\nUsage:\n"
"  trmap [-S] [-o <outfile>] [--strict-match] <ref_gff> <query_gff>\n"
"Positional arguments:\n"
"  <ref_gff>    reference annotation file name (GFF/BED format)\n"
"  <query_gff>  query file name (GFF/BED format) or \"-\" for stdin\n"
"Options:\n"
"  -o <outfile> write output to <outfile> instead of stdout\n"
"  -S           report only simple exon overlap percentages with reference\n"
"               transcripts, without classification (one line per query)\n"
"  --show-cds   add CDS:start:end info to all output transcripts\n"
"  --strict-match : when intron chains match, the '=' overlap code is assigned\n"
"               when all exons also match, otherwise assign the '~' code\n";

int main(int argc, char* argv[]) {
	GArgs args(argc, argv, "help;strict-match;show-cds;hSo:");
	args.printError(USAGE, true);
	if (args.getOpt('h') || args.getOpt("help")) {
		GMessage(USAGE);
		exit(EXIT_SUCCESS);
	}
	if (args.getOpt('S')) simpleOvl=true;
	if (args.getOpt("strict-match")) stricterMatching=true;
	if (args.getOpt("show-cds")) showCDS=true;

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
		if (t->exons.Count()==0) continue; //skip exonless entities (e.g. genes)
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
		if (!map_trees.hasKey(gseq))
			continue; //reference sequence not present in annotation, so we can't compare
		if (t->exons.Count()==0)
			continue; //only work with properly defined transcripts
		GVec<int> sidx;
		int v=0;
		sidx.Add(v); //always search the '.' strand
		if (t->strand=='+') { v=1; sidx.Add(v); }
		else if (t->strand=='-') { v=2; sidx.Add(v); }
		else { v=1; sidx.Add(v); v=2; sidx.Add(v); }
		for (int k=0;k<sidx.Count();++k) {
			GVec<GSeg*> *enu = map_trees[gseq]->it[sidx[k]].Enumerate(t->start, t->end);
			if(enu->Count()>0) {
				if (simpleOvl) {
					bool qprinted=false;
					for (int i=0; i<enu->Count(); ++i) {
						GffObj* r=(GffObj*)enu->Get(i);
						int ovlen=t->exonOverlapLen(*r);
						if (ovlen!=0) {
							float ovlcov=(100.00*ovlen)/r->len();
							if (!qprinted) {
								fprintf(outFH, "%s\t%s:%d-%d|%c", t->getID(), gseq, t->start, t->end, t->strand);
								qprinted=true;
							}
							fprintf(outFH, "\t%s:%.1f", r->getID(), ovlcov);
						}
					}
					if (qprinted) fprintf(outFH, "\n");
				} else {
					fprintf(outFH, ">%s %s:%d-%d %c ", t->getID(), t->getGSeqName(), t->start, t->end, t->strand);
					t->printExonList(outFH);
					if (showCDS && t->hasCDS()) {
					  fprintf(outFH, " CDS:");
					  t->printCDSList(outFH);
					}
					fprintf(outFH, "\n");
					for (int i=0; i<enu->Count(); ++i) {
						//static_cast<ObjInterval*>((*enu)[i])->obj->printGxf(oFile2);
						GffObj* r=(GffObj*)((*enu)[i]);
						int ovlen=0;
						char ovlcode=getOvlCode(*t, *r, ovlen, stricterMatching);
						fprintf(outFH, "%c\t", ovlcode);
						fprintf(outFH, "%s\t%c\t%d\t%d\t%s\t", r->getGSeqName(), r->strand,
							r->start, r->end, r->getID());
						r->printExonList(outFH);
						if (showCDS && r->hasCDS()) {
						  fprintf(outFH, "\tCDS:");
						  r->printCDSList(outFH);
						}
						fprintf(outFH, "\n");
					}
				}
			}
			delete enu;
		}
		delete t;
	}
	delete myQ;
    delete toFree;
	fclose(outFH);
	return 0;
}
