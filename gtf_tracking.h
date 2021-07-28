#ifndef GTF_TRACKING_H
#define GTF_TRACKING_H
/*
 *  gtf_tracking.h
 *
 */

#ifdef HEAPROFILE
#include "gperftools/heap-profiler.h"
#endif

#include "gff.h"
#include "GFaSeqGet.h"

extern int numQryFiles;
extern bool gtf_tracking_verbose;
extern bool gtf_tracking_largeScale;
extern bool qDupStrict;
extern bool stricterMatching;
extern int terminalMatchRange;
extern bool noMergeCloseExons;
extern bool debug;
extern bool reduceRefs;
//many input files, no accuracy stats are generated, no *.tmap
// and exon attributes are discarded

int cmpByPtr(const pointer p1, const pointer p2);

uint tMaxOverhang(GffObj& a, GffObj& b); //for two overlapping transcripts, return maximum terminal distance

int tMatchScore(int ovlen, GffObj* a, GffObj* b);

bool t_contains(GffObj& a, GffObj& b, bool keepAltTSS, bool intron_poking);
//returns true only IF b has fewer exons than a AND a "contains" b

char* getGSeqName(int gseq_id);

//genomic fasta sequence handling
class GFastaHandler {
 public:
  char* fastaPath;
  GFastaIndex* faIdx;
  char* getFastaFile(int gseq_id) {
     if (fastaPath==NULL) return NULL;
     Gcstr s(fastaPath);
     s.chomp('/');
     s+='/';s+=getGSeqName(gseq_id);
     Gcstr sbase(s);
     if (!fileExists(s.chars())) s.append(".fa");
     if (!fileExists(s.chars())) s.append("sta");
     if (fileExists(s.chars())) return Gstrdup(s.chars());
         else {
             GMessage("Warning: cannot find genomic sequence file %s{.fa,.fasta}\n",sbase.chars());
             return NULL;
             }
     }

   GFastaHandler(const char* fpath=NULL) {
     fastaPath=NULL;
     faIdx=NULL;
     if (fpath!=NULL && fpath[0]!=0) init(fpath);
     }

   void init(const char* fpath) {
     if (fpath==NULL || fpath[0]==0) return;
     if (!fileExists(fpath))
       GError("Error: file/directory %s does not exist!\n",fpath);
     fastaPath=Gstrdup(fpath);
     if (fastaPath!=NULL) {
         if (fileExists(fastaPath)>1) { //exists and it's not a directory
            Gcstr fainame(fastaPath);
            //the .fai name might have been given directly
            if (fainame.rindex(".fai")==fainame.length()-4) {
               //.fai index file given directly
               fastaPath[fainame.length()-4]=0;
               if (!fileExists(fastaPath))
                  GError("Error: cannot find fasta file for index %s !\n", fastaPath);
               }
              else fainame.append(".fai");
            //fainame.append(".fai");
            faIdx=new GFastaIndex(fastaPath,fainame.chars());
            Gcstr fainamecwd(fainame);
            int ip=-1;
            if ((ip=fainamecwd.rindex('/'))>=0)
               fainamecwd.cut(0,ip+1);
            if (!faIdx->hasIndex()) { //could not load index
               //try current directory
                  if (fainame!=fainamecwd) {
                    if (fileExists(fainamecwd.chars())>1) {
                       faIdx->loadIndex(fainamecwd.chars());
                       }
                    }
                  } //tried to load index
            if (!faIdx->hasIndex()) {
                 GMessage("No fasta index found for %s. Rebuilding, please wait..\n",fastaPath);
                 faIdx->buildIndex();
                 if (faIdx->getCount()==0) GError("Error: no fasta records found!\n");
                 GMessage("Fasta index rebuilt.\n");
                 FILE* fcreate=fopen(fainame.chars(), "w");
                 if (fcreate==NULL) {
                   GMessage("Warning: cannot create fasta index %s! (permissions?)\n", fainame.chars());
                   if (fainame!=fainamecwd) fcreate=fopen(fainamecwd.chars(), "w");
                   if (fcreate==NULL)
                      GError("Error: cannot create fasta index %s!\n", fainamecwd.chars());
                   }
                 if (faIdx->storeIndex(fcreate)<faIdx->getCount())
                     GMessage("Warning: error writing the index file!\n");
                 } //index created and attempted to store it
            } //multi-fasta
         } //genomic sequence given
     }
   GFaSeqGet* fetch(int gseq_id, bool checkFasta=false) {
     if (fastaPath==NULL) return NULL;
     //genomic sequence given
     GFaSeqGet* faseq=NULL;
     if (faIdx!=NULL) { //fastaPath was the multi-fasta file name
        char* gseqname=getGSeqName(gseq_id);
        GFastaRec* farec=faIdx->getRecord(gseqname);
        if (farec!=NULL) {
             faseq=new GFaSeqGet(fastaPath,farec->seqlen, farec->fpos,
                               farec->line_len, farec->line_blen);
             faseq->loadall(); //just cache the whole sequence, it's faster
             }
        else {
          GMessage("Warning: couldn't find fasta record for '%s'!\n",gseqname);
          return NULL;
          }
        }
     else //if (fileExists(fastaPath)==1)
        {
         char* sfile=getFastaFile(gseq_id);
         if (sfile!=NULL) {
            //if (gtf_tracking_verbose)
            //   GMessage("Processing sequence from fasta file '%s'\n",sfile);
            faseq=new GFaSeqGet(sfile,checkFasta);
            faseq->loadall();
            GFREE(sfile);
            }
         } //one fasta file per contig
       return faseq;
     }

   ~GFastaHandler() {
     GFREE(fastaPath);
     delete faIdx;
     }
};

bool closerRef(GffObj* a, GffObj* b, int numexons, byte rank); //for better CovLink reference ranking

class GLocus;

class COvLink {
public:
    char code;
    byte rank;
    int16_t numExons; //number of exons in the query mRNA
    GffObj* mrna;
    int ovlen;
    int16_t numJmatch; //number of matching junctions in this overlap
    COvLink(char c=0, GffObj* r=NULL, int exonCount=0, int olen=0, int jmatch=0):code(c),
    		mrna(r), ovlen(olen), numJmatch(jmatch) {
    	if (exonCount>255) exonCount=255;
    	numExons=exonCount;
		rank=classcode_rank(c);
	}

    bool operator<(COvLink& b) { //lower = higher priority
		if (rank==b.rank && b.mrna!=NULL && mrna!=NULL) {
			if (numExons==1 && b.mrna->exons.Count()!=mrna->exons.Count()) {
				if (mrna->exons.Count()==1) return true; //SET match always has priority
				else if (b.mrna->exons.Count()==1) return false;
			}
			if (numExons>1 && b.numJmatch!=numJmatch) {
				return (numJmatch > b.numJmatch);
			}
			return (ovlen==b.ovlen)? closerRef(mrna, b.mrna, numExons, rank) : (ovlen>b.ovlen);
		}
		else return rank<b.rank;
	}
    bool operator==(COvLink& b) {
		return (rank==b.rank && mrna==b.mrna);
	}
};

class GISeg: public GSeg {
 public:
   GffObj* t; //pointer to the largest transcript with a segment this exact exon coordinates
   GISeg(uint s=0,uint e=0, GffObj* ot=NULL):GSeg(s,e) { t=ot; }
};

class GIArray:public GArray<GISeg> {
  public:
   GIArray(bool uniq=true):GArray<GISeg>(true,uniq) { }
   int IAdd(GISeg* item) {
     if (item==NULL) return -1;
     int result=-1;
     if (Found(*item, result)) {
         if (fUnique) {
           //cannot add a duplicate, return index of existing item
           if (item->t!=NULL && fArray[result].t!=NULL &&
                  item->t->covlen>fArray[result].t->covlen)
               fArray[result].t=item->t;
           return result;
           }
         }
     //Found sets result to the position where the item should be
     idxInsert(result, *item);
     return result;
     }

};

class CEqList: public GList<GffObj> {
  public:
    GffObj* head;
    CEqList():GList<GffObj>((GCompareProc*)cmpByPtr, (GFreeProc*)NULL, true) {
      head=NULL;
    }
};

class CTData { //transcript associated data
public:
	GffObj* mrna; //owner transcript
	GLocus* locus;
	GList<COvLink> ovls; //overlaps with other transcripts (ref vs query)
	//GffObj* dup_of; //redundant transfrag superseded by dup_of (same query file, same locus)
	//-- just for ichain match tracking:
	GffObj* eqref; //ref transcript matching this transcript
	int qset; //qry set index (qfidx), -1 means reference dataset
	//GffObj* eqnext; //next GffObj in the linked list of matching transfrags
	bool eqhead;
	CEqList* eqlist; //keep track of matching transfrags
	//int eqdata; // flags for EQ list (is it a list head?)
	char classcode; //the best/final classcode
	// Stringtie specific data:
	double FPKM;
	double TPM;
	double cov;
	//double conf_hi;
	//double conf_lo;
	CTData(GffObj* m=NULL, GLocus* l=NULL):mrna(m), locus(l), ovls(true,true,true),
			    eqref(NULL), qset(-2), eqhead(false), eqlist(NULL),
				classcode(0), FPKM(0), TPM(0), cov(0) {
		if (mrna!=NULL) mrna->uptr=this;
	}

	~CTData() {
		ovls.Clear();
		//if ((eqdata & EQHEAD_TAG)!=0) delete eqlist;
		//if (isEqHead()) delete eqlist;
		if (eqhead) delete eqlist;
	}

  //inline bool eqHead() { return ((eqdata & EQHEAD_TAG)!=0); }
 /*  bool isEqHead() {
      if (eqlist==NULL) return false;
      return (eqlist->head==this->mrna);
      }
  */
  void joinEqList(GffObj* m) { //add list from m
   //list head is set to the transfrag with the lower qset#
  CTData* md=(CTData*)(m->uptr);
  //ASSERT(md);
  if (eqlist==NULL) { //no eqlist yet for this node
     if (md->eqlist!=NULL) { //m in an eqlist already
          eqlist=md->eqlist;
          eqlist->Add(this->mrna);
          CTData* md_head_d=(CTData*)(md->eqlist->head->uptr);
          if (this->qset < md_head_d->qset) {
               eqlist->head=this->mrna;
               eqhead=true;
               md_head_d->eqhead=false;
               }
        }
        else { //m was not in an EQ list either
          eqlist=new CEqList();
          eqlist->Add(this->mrna);
          eqlist->Add(m);
          md->eqlist=eqlist;
          if (qset<md->qset) {
        	eqlist->head=this->mrna;
        	eqhead=true;
          }
          else  {
        	eqlist->head=m;
        	md->eqhead=true;
          }
        }
      }//no eqlist before
     else { //merge two eqlists
      if (eqlist==md->eqlist) //already in the same eqlist, nothing to do
         return;
      if (md->eqlist!=NULL) {
        //copy the smaller list into the larger one
        CEqList* srclst, *destlst;
        if (md->eqlist->Count()<eqlist->Count()) {
           srclst=md->eqlist;
           destlst=eqlist;
           }
         else {
           srclst=eqlist;
           destlst=md->eqlist;
           }
         for (int i=0;i<srclst->Count();i++) {
           destlst->Add(srclst->Get(i));
           CTData* od=(CTData*)((*srclst)[i]->uptr);
           od->eqlist=destlst;
           }
        this->eqlist=destlst;
        CTData* s_head_d=(CTData*)(srclst->head->uptr);
        CTData* d_head_d=(CTData*)(destlst->head->uptr);
        if (s_head_d->qset < d_head_d->qset ) {
             this->eqlist->head=srclst->head;
             s_head_d->eqhead=true;
             d_head_d->eqhead=false;
        }
        else {
          s_head_d->eqhead=false;
          d_head_d->eqhead=true;
        }
        delete srclst;
      }
      else { //md->eqlist==NULL
        eqlist->Add(m);
        md->eqlist=eqlist;
        CTData* head_d=(CTData*)(eqlist->head->uptr);
        if (md->qset<head_d->qset) {
          eqlist->head=m;
          md->eqhead=true;
        }
      }
    }
  }

	void addOvl(TOvlData& od, GffObj* target=NULL) {
		//ovls.AddIfNew(new COvLink(code, target, ovlen));
		ovls.AddIfNew(new COvLink(od.ovlcode, target, mrna->exons.Count(), od.ovlen, od.numJmatch));
	}

	void addOvl(char code, GffObj* target=NULL, int ovlen=0) {
		ovls.AddIfNew(new COvLink(code, target, mrna->exons.Count(), ovlen));
	}

	char getBestCode(GffObj** r=NULL, int* ovlen=NULL) {
		char best_ovlcode = (ovls.Count()>0) ? ovls[0]->code : 0 ;
		if (best_ovlcode>0) {
			if (r!=NULL) *r=ovls[0]->mrna;
			if (ovlen!=NULL) *ovlen=ovls[0]->ovlen;
		}
		else {
			if (r!=NULL) *r=NULL;
			if (ovlen!=NULL) *ovlen=0;
		}
		return best_ovlcode;
    }
	bool operator<(CTData& b) { return (mrna < b.mrna); }
	bool operator==(CTData& b) { return (mrna==b.mrna); }
};


struct CEqMatch {
	int score; //match score: overlap length - overhangs
	GffObj* t;
	CTData* tdata;
	CEqMatch(GffObj* at=NULL, int sc=0):score(sc), t(at), tdata(NULL) {
		if (at!=NULL) tdata=(CTData*)(at->uptr);
	}
	bool operator<(CEqMatch& o) {
		return (score<o.score);
	}
	bool operator==(CEqMatch& o) {
	   return (score==o.score);
	}
};

class GSuperLocus;
class GTrackLocus;
class GXLocus;

class GXSeg : public GSeg {
public:
	int flags;
	GXSeg(uint s=0, uint e=0, int f=0):GSeg(s,e),flags(f) { }
};

void gatherRefLocOvls(GffObj& m, GLocus& rloc);

bool intronChainMatch(GffObj &a, GffObj &b);

//Data structure holding a query locus data (overlapping mRNAs on the same strand)
// and also the accuracy data of all mRNAs of a query locus
// (against all reference loci overlapping the same region)
class GLocus:public GSeg {
public:
    int gseq_id; //id of underlying genomic sequence
    int qfidx; // for locus tracking
    GTrackLocus* t_ptr; //for locus tracking cluster
    GffObj* mrna_maxcov;  //transcript with maximum coverage (for main "ref" transcript)
    GffObj* mrna_maxscore; //transcript with maximum gscore (for major isoform)
    GList<GffObj> mrnas; //list of transcripts (isoforms) for this locus
	GArray<GXSeg> uexons; //list of unique exons (covered segments) in this region
	GArray<GSeg> mexons; //list of merged exons in this region
	GIArray introns;
	GList<GLocus> cmpovl; //list of overlapping qry/ref loci to compare to (while forming superloci)

	//only for reference loci --> keep track of all qry loci overlaps
	//    stored in its own cmpovl list accessed by qfidx
	GPVec< GList<GLocus> > qlocovls;
	GXLocus* xlocus; //superlocus formed by exon overlaps across all qry datasets
	// -- if genomic sequence was given:
	int spl_major; // number of GT-AG splice site consensi
	int spl_rare; // number of GC-AG, AT-AC and other rare splice site consensi
	int spl_wrong; //number of "wrong" (unrecognized) splice site consensi
	int ichains; //number of multi-exon mrnas
	int ichainTP; //number of intron chains fully matching reference introns
	//int ichainATP;
	int mrnaTP;
	//int mrnaATP;
	int v; //user flag/data
	GLocus(GffObj* mrna=NULL, int qidx=-1):mrnas(true,false,false),uexons(true,true),mexons(true,true),
	  introns(), cmpovl(true,false,true), qlocovls(true) {
		//this will NOT free mrnas!
		ichains=0;
		gseq_id=-1;
		qfidx=qidx;
		t_ptr=NULL;
		creset();
		xlocus=NULL;
		mrna_maxcov=NULL;
		mrna_maxscore=NULL;
		if (mrna!=NULL) {
			start=mrna->exons.First()->start;
			end=mrna->exons.Last()->end;;
			gseq_id=mrna->gseq_id;
			GISeg seg;
			for (int i=0;i<mrna->exons.Count();i++) {
				seg.start=mrna->exons[i]->start;
				seg.end=mrna->exons[i]->end;
				int flags=0; //terminal exon flags: 1=left end, 2=right end
				if (i==0) flags|=1; //first exon
				if (i==mrna->exons.Count()-1) flags|=2; //last exon
				GXSeg xseg(seg.start, seg.end, flags);
				uexons.Add(xseg);
				mexons.Add(seg);
				if (i>0) {
					seg.start=mrna->exons[i-1]->end+1;
					seg.end=mrna->exons[i]->start-1;
					seg.t=mrna;
					introns.Add(seg);
				}
			}
			mrnas.Add(mrna);
			if (mrna->exons.Count()>1) ichains++;
			((CTData*)(mrna->uptr))->locus=this;
			mrna_maxscore=mrna;
			mrna_maxcov=mrna;
		}
	}
	void creset() {
		spl_major=0;spl_rare=0;spl_wrong=0;
		v=0; //visited/other data
		ichainTP=0;
		//ichainATP=0;
		mrnaTP=0;
		//mrnaATP=0;
		cmpovl.Clear();
	}

	void addMerge(GLocus& locus, GffObj* lnkmrna) {
		//add all the elements of the other locus (merging)
		//-- merge mexons
		GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
		int i=0; //index of first mexons with a merge
		int j=0; //index current mrna exon
		while (i<mexons.Count() && j<locus.mexons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=locus.mexons[j].start;
			uint jend=locus.mexons[j].end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//if (mexons[i].overlap(jstart, jend)) {
			//exon overlap was found :
			ovlexons.Add(j);
			//extend mexons[i] as needed
			if (jstart<istart) mexons[i].start=jstart;
			if (jend>iend) { //mexons[i] end extend
				mexons[i].end=jend;
				//now this could overlap the next mexon(s), so we have to merge them all
				while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
					uint nextend=mexons[i+1].end;
					mexons.Delete(i+1);
					if (nextend>mexons[i].end) {
						mexons[i].end=nextend;
						break; //no need to check next mexons
					}
                } //while next mexons merge
            } // mexons[i] end extend
			//  } //exon overlap
			j++; //check the next locus.mexon
		}
		//-- add the rest of the non-overlapping mexons:
		GSeg seg;
		for (int i=0;i<locus.mexons.Count();i++) {
			seg.start=locus.mexons[i].start;
			seg.end=locus.mexons[i].end;
			if (!ovlexons.Exists(i)) mexons.Add(seg);
		}
        // -- merge uexons
        //add to uexons:
		for (int i=0;i<locus.uexons.Count();i++) {
			uexons.Add(locus.uexons[i]);
		}
		for (int i=0;i<locus.introns.Count();i++) {
			introns.IAdd(&(locus.introns[i]));
            }

		// -- add locus.mrnas
		for (int i=0;i<locus.mrnas.Count();i++) {
			((CTData*)(locus.mrnas[i]->uptr))->locus=this;
			if (locus.mrnas[i]!=lnkmrna) {
				mrnas.Add(locus.mrnas[i]);
				if (locus.mrnas[i]->exons.Count()>1) ichains++;
            }
		  }
		// -- adjust start/end as needed
		if (start>locus.start) start=locus.start;
		if (end<locus.end) end=locus.end;
		if (mrna_maxcov->covlen<locus.mrna_maxcov->covlen)
			mrna_maxcov=locus.mrna_maxcov;
		if (mrna_maxscore->gscore<locus.mrna_maxscore->gscore)
			mrna_maxscore=locus.mrna_maxscore;
     }


	bool exonOverlap(GLocus& loc) {
		//check if any mexons overlap!
		int i=0;
		int j=0;
		while (i<mexons.Count() && j<loc.mexons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=loc.mexons[j].start;
			uint jend=loc.mexons[j].end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//exon overlap found
			return true;
		}
		return false;
    }

	bool add_mRNA(GffObj* mrna) {
		if (mrnas.Count()>0 && mrna->gseq_id!=gseq_id) return false; //mrna must be on the same genomic seq
		//check for exon overlap with existing mexons
		//also update uexons and mexons accordingly, if mrna is added
		uint mrna_start=mrna->exons.First()->start;
		uint mrna_end=mrna->exons.Last()->end;
		if (mrna_start>end || start>mrna_end) return false;
		bool hasovl=false;
		int i=0; //index of first mexons with a merge
		int j=0; //index current mrna exon
		GArray<int> ovlexons(true,true); //list of mrna exon indexes overlapping mexons
		while (i<mexons.Count() && j<mrna->exons.Count()) {
			uint istart=mexons[i].start;
			uint iend=mexons[i].end;
			uint jstart=mrna->exons[j]->start;
			uint jend=mrna->exons[j]->end;
			if (iend<jstart) { i++; continue; }
			if (jend<istart) { j++; continue; }
			//exon overlap found if we're here:
			ovlexons.Add(j);
			hasovl=true;
			//extend mexons[i] as needed
			if (jstart<istart) mexons[i].start=jstart;
			if (jend>iend) { //mexon stretch up
				mexons[i].end=jend;
				//now this could overlap the next mexon(s), so we have to merge them all
				while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
					uint nextend=mexons[i+1].end;
					mexons.Delete(i+1);
					if (nextend>mexons[i].end) {
						mexons[i].end=nextend;
						break; //no need to check next mexons
					}
				} //while next mexons merge
			} //possible mexons merge

			j++; //check the next mrna exon
		}//all vs all exon check loop
		if (hasovl) {
			GSeg seg;
	         //add the rest of the non-overlapping exons,
			 // and also to uexons etc.
			for (int i=0;i<mrna->exons.Count();i++) {
				seg.start=mrna->exons[i]->start;
				seg.end=mrna->exons[i]->end;
				if (!ovlexons.Exists(i)) mexons.Add(seg);
				int xterm=0;
				if (i==0) xterm|=1;
				if (i==mrna->exons.Count()-1) xterm|=2;
				GXSeg xseg(seg.start, seg.end, xterm);
				uexons.Add(xseg);
				GISeg iseg;
				if (i>0) {
					iseg.start=mrna->exons[i-1]->end+1;
					iseg.end=mrna->exons[i]->start-1;
					iseg.t=mrna;
					introns.Add(iseg);
				}
			}

			mrnas_add(mrna);
			// add to mrnas
			((CTData*)mrna->uptr)->locus=this;
			gseq_id=mrna->gseq_id;
			if (mrna->exons.Count()>1) ichains++;
		}
		return hasovl;
	}

	//simpler,basic adding of a mrna
	void mrnas_add(GffObj* mrna) {
		mrnas.Add(mrna);
		// adjust start/end
		if (start>mrna->start) start=mrna->start;
		if (end<mrna->end) end=mrna->end;
		if (mrna_maxcov->covlen<mrna->covlen) mrna_maxcov=mrna;
		if (mrna_maxscore->gscore<mrna->gscore) mrna_maxscore=mrna;
    }
};

class GSuperLocus;
class GTrackLocus;

class GSuperLocus : public GSeg {
public:
    int qfidx; //index of query dataset/file for which this superlocus was created
    GList<GLocus> qloci;
    GList<GLocus> rloci;
    GList<GffObj> qmrnas; //list of transcripts (isoforms) for this locus
    GArray<GSeg> qmexons; //list of merged exons in this region
    GArray<GXSeg> quexons; //list of unique exons (covered segments) in this region
    GIArray qintrons; //list of unique introns in this region
    //same lists for reference:
    GList<GffObj> rmrnas; //list of ref transcripts (isoforms) for this locus
    GArray<GSeg> rmexons; //list of ref merged exons in this region
    GArray<GXSeg> ruexons; //list of ref unique exons (covered segments) in this region
    GArray<GISeg> rintrons; //list of unique introns in this region
    // store problematic introns for printing:
    GIArray i_missed; //missed reference introns (not overlapped by any qry intron)
    GIArray i_notp;  //wrong ref introns (one or both ends not matching any qry intron)
    //
    GIArray i_qwrong; //totally wrong qry introns (not overlapped by any ref intron)
    GIArray i_qnotp;  //imperfect qry introns (may overlap but has no "perfect" match)


    long qbases_all;
    long rbases_all; //in fact, it's all ref bases overlapping any query loci
    int in_rmrnas; //count of ALL ref mrnas and loci given for this region
    int in_rloci; //not just those overlapping qry data
    // this will keep track of total qry loci, mrnas and exons in an area
    int total_superloci;
    int total_qloci;
    int total_qloci_alt; //total qloci with multiple transcripts

    int total_qmrnas;
    int total_qexons; //unique exons
    int total_qmexons;
    int total_qintrons; //unique introns
    int total_qichains; //total multi-exon transfrags predicted (incl. duplicates if -G)

    // NOTE: if reduceRefs these ref totals are limited to data
    //       from loci overlapping any qry loci
    int total_rmexons;
    int total_rloci;
    int total_rmrnas;
    int total_richains; //total multi-exon reference transcripts
    int total_rexons;
    int total_rintrons; //unique introns

    //--- accuracy data after compared to ref loci:
  int locusQTP;
  int locusTP; // +1 if ichainTP+mrnaTP > 0
  //int locusAQTP;
	//int locusATP; // 1 if ichainATP + mrnaATP > 0
	int locusFP;
	//int locusAFP;
	//int locusAFN;
	int locusFN;
	//---transcript level accuracy -- all exon coordinates should match (most stringent)
	int mrnaTP; // number of qry mRNAs with perfect match with ref transcripts
	//int mrnaATP;
	//---intron level accuracy (comparing the ordered set of splice sites):
	int ichainTP; // number of fully matched ref intron chains (# correctly predicted ichains)

	//int ichainFP; // number of qry intron chains not matching a reference intron chain
	//int ichainFN; // number of ref intron chains in this region not being covered by a reference intron chain
	/*
	// same as above, but Approximate -- allowing a 5bp distance around splice site coordinates
	int ichainATP; //as opposed to ichainTP, this also includes ref intron chains which are
                   //sub-chains of qry intron chains (rare cases)
     */
	//---projected features ---
	//---exon level accuracy:
	int exonTP;  //number of matched reference exons (true positives)
	int exonQTP; //number of query exons matching reference exons
	//int exonFP; //number of exons of query with no perfect match with a reference exon
	//int exonFN; //number of exons of reference with no perfect match with a query exon
	// same as the above but with acceptable approximation (10bp error window):
	/*int exonATP;
	int exonAFP;
	int exonAFN;*/

	int intronTP;  //number of perfectly overlapping introns (true positives)
	int intronFP; //number of introns of query with no perfect match with a reference intron
	int intronFN; //number of introns of reference with no perfect match with a query intron
	/*
	// same as the above but with acceptable approximation (10bp error window):
	int intronATP;
	int intronAFP;
	int intronAFN;
	*/
	//-- EGASP added these too:
	int m_exons; //number of exons totally missed (not overlapped *at all* by any query exon)
	int w_exons; //numer of totally wrong exons (query exons not overlapping *at all* any reference exon)
	int m_introns; //number of introns totally missed (not overlapped *at all* by any query intron)
	int w_introns; //numer of totally wrong introns (query introns not overlapping *at all* any reference intron)
	int m_loci; //missed loci
	int w_loci; //novel/wrong loci
	//---base level accuracy
	long baseTP; //number of overlapping bases
	long baseFP; //number of qry bases not overlapping reference
	long baseFN; //number of ref bases not overlapping qry
	//            sorted,free,unique       sorted,unique
    GSuperLocus(uint lstart=0,uint lend=0):qloci(true,false,false),rloci(true,false,false),
	qmrnas(true,false,false), qmexons(true,false), quexons(true,false), qintrons(false),
	rmrnas(true,false,false), rmexons(true,false), ruexons(true,false), rintrons(false),
	i_missed(false),i_notp(false), i_qwrong(false), i_qnotp(false){
		qfidx=-1;
		start=lstart;
		end=lend;
		qbases_all=0;
		rbases_all=0;
		baseTP=0;baseFP=0;baseFN=0;
		locusTP=0;locusQTP=0; //locusAQTP=0; locusATP=0;
		locusFP=0;// locusAFP=0;locusAFN=0;
		locusFN=0;
		in_rmrnas=0;
		in_rloci=0;
		w_loci=0;
		m_loci=0;
		total_superloci=0;
		mrnaTP=0;//mrnaFP=0;mrnaFN=0;
		ichainTP=0;//ichainFP=0;ichainFN=0;
		exonTP=0;exonQTP=0;
		//exonFP=0;exonFN=0;
		intronTP=0;intronFP=0;intronFN=0;
		/* mrnaATP=0;//mrnaAFP=0;mrnaAFN=0;
		ichainATP=0;//ichainAFP=0;ichainAFN=0;
		exonATP=0;exonAFP=0;exonAFN=0;
		intronATP=0;intronAFP=0;intronAFN=0; */
		total_rmexons=0;
		total_qmexons=0;
		total_qexons=0;total_qloci=0;total_qmrnas=0;
		total_qloci_alt=0;
		total_qintrons=0;total_qichains=0;
		total_rexons=0;total_rloci=0;total_rmrnas=0;
		total_rintrons=0;total_richains=0;
		w_exons=0;
		m_exons=0;
		w_introns=0;
		m_introns=0;
	}
    void addQlocus(GLocus& loc) {
		if (start==0 || start>loc.start) start=loc.start;
		if (end<loc.end) end=loc.end;
		qloci.Add(&loc);
		total_qloci++;
		if (loc.ichains>0 && loc.mrnas.Count()>1)
		    total_qloci_alt++;
		qmrnas.Add(loc.mrnas);
		total_qmrnas+=loc.mrnas.Count();
		total_qichains+=loc.ichains;
		qmexons.Add(loc.mexons);
		total_qmexons+=loc.mexons.Count();
		quexons.Add(loc.uexons);
		total_qexons+=loc.uexons.Count();
		qintrons.Add(loc.introns);
		total_qintrons+=loc.introns.Count();
	}

    void addRlocus(GLocus& loc) {
    	//--this is only called for ref loci that have overlaps
    	//  with at least one query locus
		if (start==0 || start>loc.start) start=loc.start;
		if (end<loc.end) end=loc.end;
		rloci.Add(&loc);
		rmrnas.Add(loc.mrnas);
		//if (reduceRefs) {
		  //partial refs counting (only for overlapping loci)
  		  total_rloci++;
		  total_rmrnas+=loc.mrnas.Count();
		  total_richains+=loc.ichains;
		  total_rmexons+=loc.mexons.Count();
		  total_rexons+=loc.uexons.Count();
		  total_rintrons+=loc.introns.Count();
		//}
		rmexons.Add(loc.mexons);
		ruexons.Add(loc.uexons);
		rintrons.Add(loc.introns);
	}

    void calcF() {
		// base level
		baseFP=qbases_all-baseTP;
		baseFN=rbases_all-baseTP;
		//exon level:
		//exonFP=total_qexons-exonTP;
		//exonFN=total_rexons-exonTP;
		//intron stats
		intronFP=total_qintrons-intronTP;
		intronFN=total_rintrons-intronTP;
		// locus/gene level:
		locusFP=total_qloci-locusQTP;
		locusFN=total_rloci-locusTP;
	}

    void addStats(GSuperLocus& s) {
		in_rmrnas+=s.in_rmrnas;
		in_rloci+=s.in_rloci;
		baseTP+=s.baseTP;
		exonTP+=s.exonTP;
		exonQTP+=s.exonQTP;
		intronTP+=s.intronTP;
		ichainTP+=s.ichainTP;
		mrnaTP+=s.mrnaTP;
		locusTP+=s.locusTP;
		locusQTP+=s.locusQTP;
		m_exons+=s.m_exons;
		w_exons+=s.w_exons;
		m_introns+=s.m_introns;
		w_introns+=s.w_introns;
		if (s.total_superloci==0 && s.qloci.Count()>0) s.total_superloci=1;
		total_superloci+=s.total_superloci;
		qbases_all+=s.qbases_all;
		rbases_all+=s.rbases_all;
		m_loci+=s.m_loci;
		w_loci+=s.w_loci;
		total_qexons+=s.total_qexons;
		total_qintrons+=s.total_qintrons;
		total_qmexons+=s.total_qmexons;
		total_rexons+=s.total_rexons;
		total_rintrons+=s.total_rintrons;
		total_rmexons+=s.total_rmexons;
		total_qmrnas+=s.total_qmrnas;
		total_qichains+=s.total_qichains;
		total_rmrnas+=s.total_rmrnas;
		total_richains+=s.total_richains;
		total_qloci+=s.total_qloci;
		total_qloci_alt+=s.total_qloci_alt;
		total_rloci+=s.total_rloci;
    }
};

class GSeqData {
	int gseq_id;
public:
    const char* gseq_name;
    GList<GffObj> refs_f; //forward strand mRNAs
    GList<GffObj> refs_r; //reverse strand mRNAs
	GList<GffObj> mrnas_f; //forward strand mRNAs
	GList<GffObj> mrnas_r; //reverse strand mRNAs
	GList<GLocus> loci_f; //forward strand loci
	GList<GLocus> loci_r; //reverse strand loci
	//--> the fields below are not used by reference data --
	GList<GSuperLocus> gstats_f; //stats for forward strand superloci
	GList<GSuperLocus> gstats_r; //stats for reverse strand superloci
	GList<GLocus> nloci_f; //"novel" loci on forward strand (no ref overlap)
	GList<GLocus> nloci_r; //"novel" loci on reverse strand (no ref overlap)
	GList<GffObj> umrnas; //unknown orientation mrnas
	GList<GLocus> nloci_u; //"novel" loci with no orientation found

	GList<CTData> tdata; //transcript data (uptr holder for all mrnas here)

	int get_gseqid() { return gseq_id; }

	//--<
	GSeqData(int gid=-1):mrnas_f(true,true,false),mrnas_r(true,true,false),
	loci_f(true,true,true),loci_r(true,true,true),
	gstats_f(true,true,false),gstats_r(true,true,false),
	nloci_f(true,false,true), nloci_r(true,false,true),
	umrnas(true,true,false), nloci_u(true,true,true), tdata(false,true,false) {
		gseq_id=gid;
		if (gseq_id>=0)
		  gseq_name=GffObj::names->gseqs.getName(gseq_id);
	}
	bool operator==(GSeqData& d){
		return (gseq_id==d.gseq_id);
	}
	bool operator>(GSeqData& d){
		return (gseq_id>d.gseq_id);
	}
	bool operator<(GSeqData& d){
		return (gseq_id<d.gseq_id);
	}
};


// a group of qry loci and a transcript cluster for a single qry dataset
class GQCluster : public GList<GffObj> {
 public:
   GffObj* mrna_maxcov;  //transcript with maximum coverage (for largest transcript)
   GffObj* mrna_maxscore; //transcript with maximum gscore ( = major isoform for Cufflinks)
   uint start;
   uint end;
   GList<GLocus> qloci;
   //GCluster cl; //just a more compact way of keeping all transcripts in these loci
   GQCluster(GList<GLocus>* loci=NULL):GList<GffObj>(true,false,false),
                                    qloci(true,false,false) {
     mrna_maxcov=NULL;
     mrna_maxscore=NULL;
     start=0;
     end=0;
     if (loci!=NULL)  {
          qloci.Add(*loci);
          for (int i=0;i<loci->Count();i++) {
             addLocus(loci->Get(i),false);
             }
          }
      }
   void addLocus(GLocus* loc, bool toLoci=true) {
     //check so we don't add locus duplicates
     if (toLoci) {
        for (int i=0;i<qloci.Count();i++) {
           if (loc==qloci[i]) return;
           }
        qloci.Add(loc);
        }
     for (int m=0;m<loc->mrnas.Count();m++) {
        GffObj* mrna=loc->mrnas[m];
        Add(mrna);
        if (start==0 || start>mrna->start) start=mrna->start;
        if (end<mrna->end) end=mrna->end;
        if (mrna_maxcov==NULL || mrna_maxcov->covlen<mrna->covlen) mrna_maxcov=mrna;
        if (mrna_maxscore==NULL || mrna_maxscore->gscore<mrna->gscore) mrna_maxscore=mrna;
        }
     }
};

//track a set of clustered qloci across multiple qry datasets
// the qloci in qcls[] overlap but not necessarily at exon level
// (so there can be multiple genes here in fact)
class GTrackLocus:public GSeg {
  public:
    char strand;
    bool hasQloci;
    //GLocus* rloc; //corresponding reference locus, if available
    GList<GLocus> rloci; //ref loci found overlapping this region
    GVec<GQCluster*> qcls; //all qloci for this superlocus, grouped by dataset
    GTrackLocus(int numqryfiles, GLocus* qloc=NULL, int q=-1):GSeg(0,0),rloci(true,false,true),qcls() {
      strand='.';
      if (numqryfiles>0) {
    	  qcls.Resize(numqryfiles, NULL);
      }
      else GError("Error: invalid GTrackLocus constructor called before initializing numQueryFiles.\n");
      if (qloc!=NULL) addQLocus(qloc,q);
      }
    void init(int num) {
        //for (int i=0;i<num;i++) qcls.Add(NULL);
    	qcls.Resize(num, NULL);
    }

    void addRLocus(GLocus* rl) {
      if (rl==NULL) return;
      if (rl->qfidx>=0)
          GError("Error: GTrackLocus::addRLocus called with a query locus (set# %d)\n",
                        rl->qfidx+1);
      if (strand=='.') strand=rl->mrna_maxcov->strand;
      if (start==0 || start>rl->start) start=rl->start;
      if (end==0 || end<rl->end) end=rl->end;
      rl->t_ptr=this;
      rloci.Add(rl);
      }

    void addQLocus(GLocus* loc, int q=-1) { //adding qry locus
      if (loc==NULL) return;
      if (strand=='.' && loc->mrna_maxcov->strand!='.')
           strand=loc->mrna_maxcov->strand;
      if (loc->qfidx<0 && q<0)
         GError("Error at GTrackLocus::addQLocus(): locus.qfidx not set and index not given!\n");
      if (q>=0) loc->qfidx=q;
           else q=loc->qfidx;
      if (start==0 || start>loc->start) start=loc->start;
      if (end==0 || end<loc->end) end=loc->end;
      if (qcls[q]==NULL) qcls[q]=new GQCluster();
      hasQloci=true;
      loc->t_ptr = this;
      qcls[q]->addLocus(loc);
      }

    bool add_Locus(GLocus* loc) {
      if (start==0 || overlap(*loc)) { //simple range overlap, not exon overlap
         if (loc->qfidx<0) addRLocus(loc);
                      else addQLocus(loc);
         return true;
         }
      return false;
      }


   void addQCl(int q, GQCluster* qcl, GLocus* lnkloc) {
      for (int i=0;i<qcl->qloci.Count();i++) {
         GLocus* loc=qcl->qloci[i];
         if (loc==lnkloc) continue; // or if loc->t_ptr==this ?
         hasQloci=true;
         loc->t_ptr=this;
         qcls[q]->addLocus(loc);
         }
     }

   void addMerge(GTrackLocus* loctrack, int qcount, GLocus* lnkloc) {
      if (loctrack==NULL) return;
      //merge qloci
      for (int q=0; q < qcount; q++) {
          if (qcls[q]==NULL) {
             if (loctrack->qcls[q]!=NULL) {
                 qcls[q]=loctrack->qcls[q];
                 loctrack->qcls[q]=NULL; //just move pointer here
                 //set all t_ptr pointers for moved loci
                 for (int ql = 0; ql < qcls[q]->qloci.Count(); ql++) {
                    qcls[q]->qloci[ql]->t_ptr=this;
                    }
                 hasQloci=true;
                 }
             }
           else //existing qloci at q
             if (loctrack->qcls[q]!=NULL) { //merge elements
                addQCl(q, loctrack->qcls[q], lnkloc);
              }
          }//for each qset
      //merge rloci, if any
      if (loctrack->rloci.Count()>0) {
         for (int l=0;l<loctrack->rloci.Count();l++) {
           if (loctrack->rloci[l]!=lnkloc && loctrack->rloci[l]->t_ptr!=this) {
              rloci.Add(loctrack->rloci[l]);
              loctrack->rloci[l]->t_ptr=this;
              }
           }
         }
      if (loctrack->start<start) start=loctrack->start;
      if (loctrack->end>end) end=loctrack->end;
      if (strand=='.' && loctrack->strand!='.') strand=loctrack->strand;
      }

    /*
    void add_QLoci(GList<GLocus>* loci, int q, GLocus& r) {
      // only add loci overlapping given refloc
      //rloc=&r;
      if (loci==NULL) return;
      for (int i=0;i<loci->Count();i++) {
         GLocus* loc=loci->Get(i);
         // if (!loc->exonOverlap(r)) continue;  //do we really needed exon overlap?
         if (!loc->overlap(r)) continue;
         if (start==0 || start>loc->start) start=loc->start;
         if (end==0 || end<loc->end) end=loc->end;
         loc->t_ptr=this;
         loc->qfidx=q;
         if (qcls[q]==NULL) qcls[q]=new GQCluster();
         qcls[q]->addLocus(loc);
         }
      strand=r.mrnas[0]->strand;
      }
     */
    ~GTrackLocus() {
      for (int q=0;q<numQryFiles;q++)
           if (qcls[q]!=NULL) { delete qcls[q]; qcls[q]=NULL; }
      }

    GQCluster* operator[](int q) {
      if (q<0 || q>=numQryFiles)
          GError("Error: qfidx index out of bounds (%d) for GTrackLocus!\n",q);
      return qcls[q];
      }
};

class GXConsensus:public GSeg {
 public:
   static int count;
   int id; //XConsensus ID
   int tss_id; //group id for those xconsensi with shared first exon
   int p_id; //group id for those xconsensi with "similar" protein
   GffObj* tcons; //longest transcript to represent the combined "consensus" structure
   GffObj* ref; //overlapping reference transcript
   char refcode; // the code for ref relationship (like in the tracking file)
   char* aa;
   int aalen;
   GXConsensus* contained; //if contained into another GXConsensus
   //list of ichain-matching query transcripts that contributed to this consensus
   GList<GffObj> qchain;
   GXConsensus(GffObj* c, CEqList* qlst, GffObj* r=NULL, char rcode=0)
                   :qchain(false,false,false) {
      ref=r;
      refcode=rcode;
      tcons=c;
      if (qlst!=NULL) qchain.Add(*((GList<GffObj>*)qlst));
                 else qchain.Add(c);
      count++;
      tss_id=0;
      p_id=0;
      aalen=0;
      id=count;
      aa=NULL;
      start=tcons->start;
      end=tcons->end;
      contained=NULL;
      }
   ~GXConsensus() {
     if (aa!=NULL) GFREE(aa);
     }
};

//cross sample and reference superlocus data structure
class GXLocus:public GSeg {
 public:
    int id;
    int num_mtcons; //number of multi-exon "consensus" transcripts in this locus
    char strand;
    GList<GLocus> rloci; //list of ref loci overlapping any of the mexons
    GList<GLocus> qloci; //loci from all qry datasets that have overlapping exons with this region
    GArray<GSeg> mexons; //list of merged exonic regions for this locus
    GList<GXConsensus> tcons;
    GXLocus(GLocus* lfirst=NULL):GSeg(0,0),
        rloci((GCompareProc*)cmpByPtr, (GFreeProc*)NULL, true),
        qloci((GCompareProc*)cmpByPtr, (GFreeProc*)NULL, true),
        mexons(true,true), tcons(true,true,false) {
      strand='.';
      num_mtcons=0;
      if (lfirst!=NULL) {
         add_Locus(lfirst);
         }
      id=0;
      }

    bool add_Locus(GLocus* loc) {
      if (mexons.Count()>0 && (end<loc->start || start > loc->end))
              return false; //no chance for overlapping exons
      if (mexons.Count()==0) {
          mexons.Add(loc->mexons);
          start=loc->start;
          end=loc->end;
          if (loc->qfidx<0) rloci.Add(loc);
                      else  qloci.Add(loc);
          strand=loc->mrna_maxcov->strand;
          loc->xlocus=this;
          return true;
          }
      int f=0;
      if (loc->qfidx<0) {
        if (rloci.Found(loc,f)) return false;
        }
      else if (qloci.Found(loc,f)) return false;

      // -- merge mexons
      GArray<int> ovlexons(true,true); //list of locus.mexons indexes overlapping existing mexons
      int i=0; //index of first mexons with a merge
      int j=0; //index current mrna exon
      while (i<mexons.Count() && j<loc->mexons.Count()) {
          uint istart=mexons[i].start;
          uint iend=mexons[i].end;
          uint jstart=loc->mexons[j].start;
          uint jend=loc->mexons[j].end;
          if (iend<jstart) { i++; continue; }
          if (jend<istart) { j++; continue; }
          //if (mexons[i].overlap(jstart, jend)) {
          //exon overlap was found :
          ovlexons.Add(j);
          //extend mexons[i] as needed
          if (jstart<istart) mexons[i].start=jstart;
          if (jend>iend) { //mexons[i] end extend
              mexons[i].end=jend;
              //now this could overlap the next mexon(s), so we have to merge them all
              while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                  uint nextend=mexons[i+1].end;
                  mexons.Delete(i+1);
                  if (nextend>mexons[i].end) {
                      mexons[i].end=nextend;
                      break; //no need to check next mexons
                  }
              } //while next mexons merge
          } // mexons[i] end extend
          //  } //exon overlap
          j++; //check the next locus.mexon
        }//while mexons
      if (ovlexons.Count()==0) return false;
      if (strand=='.' && loc->mrna_maxcov->strand!='.')
             strand=loc->mrna_maxcov->strand;
      //have exon overlap:
      //-- add the rest of the non-overlapping mexons:
       GSeg seg;
       for (int i=0;i<loc->mexons.Count();i++) {
            seg.start=loc->mexons[i].start;
            seg.end=loc->mexons[i].end;
            if (!ovlexons.Exists(i)) mexons.Add(seg);
            }
      // -- adjust start/end as needed
      if (start>loc->start) start=loc->start;
      if (end<loc->end) end=loc->end;
      loc->xlocus=this;
      if (loc->qfidx<0) rloci.Add(loc);
                  else  qloci.Add(loc);
      return true;
      }

  void addMerge(GXLocus& oxloc) {
    GArray<int> ovlexons(true,true); //list of oxloc.mexons indexes overlapping existing mexons
    int i=0; //index of first mexons with a merge
    int j=0; //index current mrna exon
    while (i<mexons.Count() && j<oxloc.mexons.Count()) {
        uint istart=mexons[i].start;
        uint iend=mexons[i].end;
        uint jstart=oxloc.mexons[j].start;
        uint jend=oxloc.mexons[j].end;
        if (iend<jstart) { i++; continue; }
        if (jend<istart) { j++; continue; }
        //if (mexons[i].overlap(jstart, jend)) {
        //exon overlap was found :
        ovlexons.Add(j);
        //extend mexons[i] as needed
        if (jstart<istart) mexons[i].start=jstart;
        if (jend>iend) { //mexons[i] end extend
            mexons[i].end=jend;
            //now this could overlap the next mexon(s), so we have to merge them all
            while (i<mexons.Count()-1 && mexons[i].end>mexons[i+1].start) {
                uint nextend=mexons[i+1].end;
                mexons.Delete(i+1);
                if (nextend>mexons[i].end) {
                    mexons[i].end=nextend;
                    break; //no need to check next mexons
                }
            } //while next mexons merge
        } // mexons[i] end extend
        //  } //exon overlap
        j++; //check the next oxloc.mexon
    }
    if (ovlexons.Count()==0) {
      GError("Error: attempt to merge GXLoci with non-overlapping exons!\n");
      }
    //-- add the rest of the non-overlapping mexons:
    GSeg seg;
    for (int i=0;i<oxloc.mexons.Count();i++) {
        seg.start=oxloc.mexons[i].start;
        seg.end=oxloc.mexons[i].end;
        if (!ovlexons.Exists(i)) mexons.Add(seg);
        }
   if (start>oxloc.start) start=oxloc.start;
   if (end<oxloc.end) end=oxloc.end;
   if (strand=='.') strand=oxloc.strand;
    //-- steal all qloci and rloci
    for (int i=0;i<oxloc.qloci.Count();i++) {
         if (oxloc.qloci[i]->xlocus==this) continue;
         qloci.Add(oxloc.qloci[i]);
         oxloc.qloci[i]->xlocus=this;
         }
    for (int i=0;i<oxloc.rloci.Count();i++) {
         if (oxloc.rloci[i]->xlocus==this) continue;
         rloci.Add(oxloc.rloci[i]);
         oxloc.rloci[i]->xlocus=this;
         }
  } //::addMerge()


 void checkContainment(bool keepAlt5=false, bool intron_poking=false) {
   //checking containment
  for (int j=0;j<tcons.Count()-1;j++) {
    GXConsensus* t=tcons[j];
    for (int i=j+1;i<tcons.Count();i++) {
       //if (tcons[i]->contained!=NULL && t->tcons->exons.Count()>1) continue; //will check the container later anyway
       int c_status=checkXConsContain(t->tcons, tcons[i]->tcons, keepAlt5, intron_poking);
       if (c_status==0) continue; //no containment relationship between t and tcons[i]
       if (c_status>0) { //t is a container for tcons[i]
            tcons[i]->contained=t;
            }
         else { //contained into exising XCons
            t->contained=tcons[i];
            break;
            }
       }
   }
  }

 int checkXConsContain(GffObj* a, GffObj* b, bool keepAltTSS, bool intron_poking) {
  // returns  1 if a is the container of b
  //         -1 if a is contained in b
  //          0 if no
  if (a->end<b->start || b->end<a->start) return 0;
  if (a->exons.Count()==b->exons.Count()) {
      if (a->exons.Count()>1) {
    	  /*if (((CTData*)a->uptr)->qset!=((CTData*)b->uptr)->qset)
    	     return 0; //different sample, same number of exons - no containment possible
    	               //because equivalence was already tested across samples

    	  */
    	  if (intronChainMatch(*a, *b)) {
        	 //if matched choose one to be the "container"
        	 //TODO: ideally we should choose one
        	 //    a) the one better represented across multiple samples
        	 //       OR
        	 //    b) the one with a closer match of terminal exons to a reference transcript
        	 // for now we choose the longest one (which might not be the "right one")
    		 return ((a->covlen > b->covlen) ? 1 : -1) ;
    		 }
    	  }
      else { //single exon containment testing
             //this is fuzzy and messy (end result may vary depending on the testing order)
             int ovlen=a->exons[0]->overlapLen(b->exons[0]);
             int lmax=a->covlen;
             int lmin=b->covlen;
             if (lmin>lmax) Gswap(lmin,lmax);
             //if (ovlen>=lmax*0.7 || ovlen>=lmin*0.85) { //if at least 85% of the shorter one is covered, it is contained
             if (ovlen>=lmin*0.85) { //if at least 85% of the shorter one is covered, it is considered contained
                return ((a->covlen>b->covlen) ? 1 : -1);
                }
              else return 0;
           }
     }
   //different number of exons:
   if (a->exons.Count()>b->exons.Count())
	    return t_contains(*a, *b, keepAltTSS, intron_poking) ?  1 : 0;
   else return t_contains(*b, *a, keepAltTSS, intron_poking) ? -1 : 0;
  }
 /*
 void addXCons(GXConsensus* t) {
  tcons.Add(t);
 }
 */

}; //GXLocus



int parse_mRNAs(GfList& mrnas,
				 GList<GSeqData>& glstdata,
				 bool is_ref_set=true,
				 bool discardDups=false,
				 int qfidx=-1, bool only_multiexon=false);

//reading a mRNAs from a gff file and grouping them into loci
void read_mRNAs(FILE* f, GList<GSeqData>& seqdata, GList<GSeqData>* ref_data=NULL,
              bool discardDups=false, int qfidx=-1, const char* fname=NULL,
              bool only_multiexon=false);

void read_transcripts(FILE* f, GList<GSeqData>& seqdata,
#ifdef CUFFLINKS
  boost::crc_32_type& crc_result,
#endif
  bool keepAttrs=true);

void sort_GSeqs_byName(GList<GSeqData>& seqdata);

//moved into gff.h:
//char singleExonTMatch(GffObj& m, GffObj& r, int& ovlen);

/*
//strict intron chain match, or single-exon match
bool tMatch(GffObj& a, GffObj& b, int& ovlen, bool relaxed_singleExonMatch=false,
           bool contain_only=false);
*/

//use qsearch to "position" a given coordinate x within a list of transcripts sorted
//by their start (lowest) coordinate;
//the return value is the index of the closest GffObj starting just *ABOVE* coordinate x
//Convention: returns -1 if there is no such GffObj (i.e. mrnas.Last().start <= x)
int qsearch_mrnas(uint x, GList<GffObj>& mrnas);
int qsearch_loci(uint x, GList<GLocus>& segs); // same as above, but searching for loci segments
//--- qsearch transcript overlap (returning -1 when no overlap is found)

GSeqData* getRefData(int gid, GList<GSeqData>& ref_data); //returns reference GSeqData for a specific genomic sequence

#endif
