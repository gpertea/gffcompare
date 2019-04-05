#include <iostream>

#include <fstream>
#include <sstream>

#include "GArgs.h"
#include "gff.h"
#include "GIntervalTree.h"

using std::cout;

#define VERSION "0.10.6"

bool simpleOvl=false;
bool strictMatching=false;

struct GSTree {
	GIntervalTree it[3]; //0=unstranded, 1: + strand, 2 : - strand
};

int main(int argc, char* argv[]) {
	const std::string usage = std::string("trmap v" VERSION)+
	" : transcript to reference mapping and overlap classifier.\nUsage:\n"+
	"  trmap [-S] [-o <outfile>] [--strict-match] <ref_gff> <query_gff>\n"+
	        "Positional arguments:\n"+
			"  <ref_gff>    reference annotation file name (GFF/BED format)\n"+
			"  <query_gff>  query file name (GFF/BED format) or \"-\" for stdin\n"+
			"Options:\n"+
			"  -o <outfile> write output to <outfile> instead of stdout\n"+
			"  -S           report only simple reference overlap percentages, without\n"+
			"               classification (one line per query)\n"+
			"  --strict-match : when intron chains match, the '=' overlap code is assigned\n"+
			"               when all exons also match, otherwise assign the '~' code\n";
	GArgs args(argc, argv, "help;strict-match;hSo:");
	args.printError(usage.c_str(), true);
	if (args.getOpt('h') || args.getOpt("help")) {
		cout << usage;
		exit(EXIT_SUCCESS);
	}
	if (args.getOpt('S')) simpleOvl=true;
    if (args.getOpt("strict-match")) strictMatching=true;
	GHash<GSTree> map_trees;

	const char* o_file = args.getOpt('o') ? args.getOpt('o') : "-";

	if (args.startNonOpt()!=2) {
		std::cerr << usage << "\nOnly " << args.startNonOpt() << " arguments provided (expected 2)\n";
		exit(1);
	}
	const char* ref_file = args.nextNonOpt();
	const char* q_file = args.nextNonOpt();

	FILE* fr=fopen(ref_file, "r");

	//always good to check if the file is actually there and can be read
	if (fr==NULL) GError("Error: could not open reference annotation file (%s)!\n", ref_file);

	GffReader myR(fr, true, true);
	const char* fext=getFileExt(ref_file);
	if (Gstricmp(fext, "bed")==0) myR.isBED();
	GffObj* t=NULL;
	GPVec<GffObj> toFree(true);
	while ((t=myR.readNext())!=NULL) {
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
		toFree.Add(t);
	}
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
	GffReader myQ(fq, true, true);
	if (fext && Gstricmp(fext, "bed")==0) myQ.isBED();
	//myQ.readAll(false, true, true);
	//for (int i=0; i<myQ.gflst.Count(); i++) {
	//	GffObj* t=myQ.gflst[i];
	t=NULL;
	while ((t=myQ.readNext())!=NULL) {
		//if (map_trees.count(t->getGSeqName())==0) continue;
		const char* gseq=t->getGSeqName();
		if (!map_trees.hasKey(gseq)) continue;
		GVec<int> sidx;
		int v=0;
		sidx.Add(v); //always search the '.' strand
		if (t->strand=='+') { v=1; sidx.Add(v); }
		else if (t->strand=='-') { v=2; sidx.Add(v); }
		else { v=1; sidx.Add(v); v=2; sidx.Add(v); }
		for (int k=0;k<sidx.Count();++k) {
			TemplateStack<GSeg*> * enu = map_trees[gseq]->it[sidx[k]].Enumerate(t->start, t->end);
			if(enu->Size()!=0) {
				if (simpleOvl) {
					fprintf(outFH, "%s\t%s:%d-%d|%c", t->getID(), gseq, t->start, t->end, t->strand);
					for (int i=0; i<enu->Size(); ++i) {
						//static_cast<ObjInterval*>((*enu)[i])->obj->printGxf(oFile2);
						GffObj* r=(GffObj*)((*enu)[i]);
						int ovlen=t->overlapLen(r);
						if (ovlen==0)
							GError("Error: zero length simple overlap reported! (%s vs %s)\n", t->getID(), r->getID());
						float ovlcov=(100.00*ovlen)/r->len();
						fprintf(outFH, "\t%s:%.1f", r->getID(), ovlcov);
						//if (i+1<enu->Size()) fprintf(outFH, ",");
					}
					fprintf(outFH, "\n");
				} else {
					fprintf(outFH, ">%s %s:%d-%d %c ", t->getID(), t->getGSeqName(), t->start, t->end, t->strand);
					t->printExonList(outFH);
					fprintf(outFH, "\n");
					for (int i=0; i<enu->Size(); ++i) {
						//static_cast<ObjInterval*>((*enu)[i])->obj->printGxf(oFile2);
						GffObj* r=(GffObj*)((*enu)[i]);
						int ovlen=0;
						char ovlcode=getOvlCode(*t, *r, ovlen, strictMatching);
						fprintf(outFH, "%c\t", ovlcode);
						r->printGTab(outFH);
					}
				}
			}
			delete enu;
		}
		delete t;
	}

	fclose(outFH);
	return 0;
}
