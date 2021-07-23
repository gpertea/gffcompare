#include "gtf_tracking.h"

bool gtf_tracking_verbose = false;
bool gtf_tracking_largeScale=false; //many input Cufflinks files processed at once by cuffcompare, discard exon attributes

int numQryFiles=0;

bool reduceRefs=false; //-R

bool qDupStrict=false;
bool stricterMatching=false;
int terminalMatchRange=0;
bool noMergeCloseExons=false;
bool debug=false;

int GXConsensus::count=0;

char* getGSeqName(int gseq_id) {
 return GffObj::names->gseqs.getName(gseq_id);
}

int cmpByPtr(const pointer p1, const pointer p2) {
  return (p1>p2) ? 1: ((p1==p2)? 0 : -1);
  }

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

//for two overlapping transcripts, return maximum terminal distance
uint tMaxOverhang(GffObj& a, GffObj& b){
	//WARNING: this does not make sense if a and b do not overlap!
	uint dstart=(a.start>b.start) ? a.start-b.start : b.start-a.start;
	uint dend=(a.end>b.end) ? a.end-b.end : b.end-a.end;
    return ((dstart>dend) ? dstart : dend);
}

int tMatchScore(int ovlen, GffObj* a, GffObj* b) { //simply ovlen - overhangs
	//WARNING: this does not make sense if a and b do not overlap!
	int dstart=(a->start>b->start) ? a->start-b->start : b->start-a->start;
	int dend=(a->end>b->end) ? a->end-b->end : b->end-a->end;
    return ((int)ovlen - dstart - dend); //can be negative for large overhangs
}

GffObj* is_TDup(GffObj* m, GList<GffObj>& mrnas, int& dupidx, bool matchContain=false) {
 //mrnas MUST be sorted by start coordinate
 //this is optimized for when mrnas list is being populated, in sorted order
 //as it starts scanning from the end of the list
  dupidx=-1;
  if (mrnas.Count()==0) return NULL;
  //int nidx=qsearch_mrnas(m->end, mrnas);
  //if (nidx==0) return NULL;
  //if (nidx==-1) nidx=mrnas.Count();//all can overlap
  int nidx=mrnas.Count();
  for (int i=nidx-1;i>=0;i--) {
      GffObj& omrna=*mrnas[i];
      if (m->start>omrna.end) {
           if (m->start-omrna.start>GFF_MAX_EXON) break; //give up already, went too far back
           continue;
           }
      if (omrna.start>m->end) continue; //this should never be the case if nidx was found with qsearch_mrnas(m->end)
      //locus overlap here:
      //if (tMatch(*m, omrna, ovlen, !matchContain, matchContain)) {
      int ovlen=0;
      char matchType=transcriptMatch(*m, omrna, ovlen);
      if (matchType>0) {
    	  if (matchType=='=' || !matchContain || omrna.contains(m) ) {
				 dupidx=i;
				 return mrnas[i];
		  }
      }
  }
  return NULL;
}

bool intronChainMatch(GffObj&a, GffObj&b) {
	if (a.exons.Count()!=b.exons.Count()) return false;
	if (a.exons.Count()<2) return false;
 	for (int i=1;i<a.exons.Count();i++) {
 		//if (i<imax) ovlen+=a.exons[i]->len();
 		if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
 				(a.exons[i]->start!=b.exons[i]->start)) {
 			return false; //intron mismatch
 		}
 	}
 	return true;
}

bool intronRedundant(GffObj& ti, GffObj&  tj, bool no5diff=false, bool intron_poking=false) {
	//two transcripts are "intron redundant" iff one transcript's intron chain
	// is a sub-chain of the other's
	//no5diff=true : will NOT consider redundant if they have a different first intron at 5'
	//allowXintron=true : allow a contained transcript to start or end within a container's intron
	int imax=ti.exons.Count()-1;
	int jmax=tj.exons.Count()-1;
	if (imax==0 || jmax==0) return false; //don't deal with single-exon transcripts here
	if (ti.exons[imax]->start<tj.exons[0]->end ||
			tj.exons[jmax]->start<ti.exons[0]->end )
		return false; //intron chains do not overlap at all

	uint eistart=0, eiend=0, ejstart=0, ejend=0; //exon boundaries
	int i=1; //exon idx to the right of the current intron of ti
	int j=1; //exon idx to the right of the current intron of tj
	//find the first intron overlap:
	while (i<=imax && j<=jmax) {
		eistart=ti.exons[i-1]->end;
		eiend=ti.exons[i]->start;
		ejstart=tj.exons[j-1]->end;
		ejend=tj.exons[j]->start;
		if (ejend<eistart) { j++; continue; }
		if (eiend<ejstart) { i++; continue; }
		//we found an intron overlap
		break;
	}
	if (eistart!=ejstart || eiend!=ejend)
		return false; //first intron overlap is NOT an exact intron match

	if ((i>1 && j>1) //not the first intron for at least one of the transcripts
			|| i>imax || j>jmax) { //no intron overlap found
		return false;
	}
	//we have the first matching intron on the left
	if (!intron_poking) {
		if (j>i //ti starts within tj (ti probably contained within tj)
				//i==1, ti's start must not conflict with the previous intron of tj
				&& ti.start<tj.exons[j-1]->start) return false;
		//or tj contained within ti?
		if (i>j && tj.start<ti.exons[i-1]->start) return false;
	}
	// ---- comment out the next 2 if statements below if just "intron compatibility"
	//      (i.e. extension of intron chains) is desired
	if (j>i  //ti starts within tj (ti probably contained within tj)
		 && // then tj must contain ti, so ti's last intron must end with or before tj's last intron
		 ti.exons[imax]->start>tj.exons[jmax]->start) return false;
	if (i>j && tj.exons[jmax]->start>ti.exons[imax]->start) return false;
	// ----

	//now check if the rest of the introns match in the same sequence
	int i_start=i; //first (leftmost) matching intron of ti (1-based index)
	int j_start=j; //first (leftmost) matching intron of tj
	i++;j++;
	while (i<=imax && j<=jmax) {
		if (ti.exons[i-1]->end!=tj.exons[j-1]->end ||
				ti.exons[i]->start!=tj.exons[j]->start) return false;
		i++;
		j++;
	}
	i--; j--; //i,j=indexes of last (rightmost) matching intron i_end, j_end
	if (!intron_poking) { //check for terminal exons of the contained sticking out within a container's intron
		if (i==imax && j<jmax && //tj has more introns
			// check if ti's end doesn't conflict with the current tj exon boundary
			 ti.end>tj.exons[j]->end) return false;

		if (j==jmax && i<imax &&
			tj.end>ti.exons[i]->end) return false;
	}
	if (no5diff && imax!=jmax) { //different number of introns
		//if they start with a different 5' intron
		// they are NOT considered redundant
		if (ti.strand=='+') {
			if (i_start!=j_start) return false;
		}
		else { //reverse strand
			if (imax-i!=jmax-j) return false;
		}
	}
	return true; //they are intron-redundant
}


bool t_contains(GffObj& a, GffObj& b, bool keepAltTSS, bool intron_poking) {
 //returns true if b's intron chain (or single exon) is included in a
 if (b.exons.Count()>=a.exons.Count()) return false;
 if (b.exons.Count()==1) {
    //check if b is contained in any of a's exons:
    for (int i=0;i<a.exons.Count();i++) {
       if (b.start>=a.exons[i]->start && b.end<=a.exons[i]->end) return true;
       }
    return false;
    }

 if (intronRedundant(a, b, keepAltTSS, intron_poking)) {
     //intronRedudant allows b's initial/terminal exons to extend beyond a's boundaries
     //but we don't allow this here *unless* user already relaxed the redundancy conditions!
     if (intron_poking) return true;
     else return (b.start>=a.start && b.end<=a.end);
    }
  else return false;
 }


int is_Redundant(GffObj*m, GList<GffObj>* mrnas, bool no5share=false, bool intron_poking=false) {
	//first locate the list index of the mrnas starting just ABOVE m->end
	if (mrnas->Count()==0) return -1;
	int nidx=qsearch_mrnas(m->end, *mrnas);
	if (nidx==0) return -1; //none can overlap
	if (nidx==-1) nidx=mrnas->Count();//all can overlap
	for (int i=nidx-1;i>=0;i--) {
		GffObj& t=*mrnas->Get(i); //overlap check target
		if (m->start>t.end) { //m starts after target ends
			if (m->start > t.start+GFF_MAX_LOCUS)
				break; //went too far back, give up
			continue;
		}
		if (t.start>m->end) continue; //this should never be the case if nidx was found correctly
		//what about single-exon transcript redundancy ? probably not needed within a sample
		if (intronRedundant(*m, t, no5share, intron_poking)) return i;
	}
	return -1;
}

bool t_dominates(GffObj* a, GffObj* b) {
 // for redundant / intron compatible transfrags:
 // returns true if a "dominates" b, i.e. a has more exons or is longer
 if (a->exons.Count()==b->exons.Count())
         return (a->covlen>b->covlen);
    else return (a->exons.Count()>b->exons.Count());
}

bool betterTDup(GffObj* a, GffObj* b) {
  if (a->exons.Count()!=b->exons.Count())
    return (a->exons.Count()>b->exons.Count());
  if (a->hasCDS()!=b->hasCDS())
     return (a->hasCDS()>b->hasCDS());
   //for annotation purposes, it's more important to keep the
   //longer transcript, instead of the one that was loaded first
  if (a->covlen != b->covlen)
         return (a->covlen > b->covlen);
    else return (a->track_id < b->track_id);
}

int parse_mRNAs(GfList& mrnas,
				 GList<GSeqData>& glstdata,
				 bool is_ref_set, bool discardDups,
				 int qfidx, bool only_multiexon) {
	int tredundant=0; //redundant transcripts discarded
	int total_kept=0;
	//int total_seen=mrnas.Count();
	for (int k=0;k<mrnas.Count();k++) {
		GffObj* m=mrnas[k];
		int i=-1;
		GSeqData f(m->gseq_id);
		GSeqData* gdata=NULL;
		uint tlen=m->len();
		if (m->hasErrors() || (tlen+500>GFF_MAX_LOCUS)) { //should probably report these in a file too..
			if (gtf_tracking_verbose)
			      GMessage("Warning: transcript %s discarded (structural errors found, length=%d).\n", m->getID(), tlen);
			continue;
			}
		if (only_multiexon && m->exons.Count()<2) {
			continue;
			}
		//GStr feature(m->getFeatureName());
		//feature.lower();
		//bool gene_or_locus=(feature.endsWith("gene") ||feature.index("loc")>=0);
		//if (m->exons.Count()==0 && gene_or_locus) {
		if (m->isDiscarded()) {
			//discard generic "gene" or "locus" features with no other detailed subfeatures
			if (!is_ref_set && gtf_tracking_verbose)
			   GMessage("Warning: discarding non-transfrag (GFF generic gene/locus container?) %s\n",m->getID());
			continue;
			}

		if (m->exons.Count()==0) {
			if (gtf_tracking_verbose && !is_ref_set)
			    GMessage("Warning: %s %s found without exon segments (adding default exon).\n",m->getFeatureName(), m->getID());
			m->addExon(m->start,m->end);
		}
		if (glstdata.Found(&f,i)) gdata=glstdata[i];
		else {
			gdata=new GSeqData(m->gseq_id);
			glstdata.Add(gdata);
		}

		double fpkm=0;
		double cov=0;
		double tpm=0;
		//GffObj* dup_by=NULL;
		GList<GffObj>* target_mrnas=NULL;
		if (is_ref_set) { //-- ref transcripts
		   if (m->strand=='.') {
		     //unknown strand - discard from reference set (!)
			 if (gtf_tracking_verbose)
				 GMessage("Warning: reference transcript %s has undetermined strand, discarded.\n", m->getID());
		     continue;
		   }
		   total_kept++;
		   target_mrnas=(m->strand=='+') ? &(gdata->mrnas_f) : &(gdata->mrnas_r);
		   if (discardDups) {
		     //check all gdata->mrnas_r (ref_data) for duplicate ref transcripts
		     int rpidx=-1;
		     GffObj* rp= is_TDup(m, *target_mrnas, rpidx, true);
		       //always strict checking of reference duplicates: containment required
		     if (rp!=NULL) { //duplicate found
		      //discard one of them
		      //but let's keep the gene_name if present
		      //DEBUG:
		      //GMessage("Ref duplicates: %s = %s\n", rp->getID(), m->getID());
		      tredundant++;
		      total_kept--;
		      if (betterTDup(rp, m)) {
		           if (rp->getGeneName()==NULL && m->getGeneName()!=NULL) {
		                  rp->setGeneName(m->getGeneName());
		           }
				  if (debug)
					   GMessage("\tReference transcript %s discarded (duplicate of %s)\n",
					      m->getID(), rp->getID() );
		           continue;
		      }
		      else {
		           if (m->getGeneName()==NULL && rp->getGeneName()!=NULL) {
		                  m->setGeneName(rp->getGeneName());
		           }
				  if (debug)
					   GMessage("\tReference transcript %s discarded (duplicate of %s)\n",
					      rp->getID(), m->getID() );
		           ((CTData*)(rp->uptr))->mrna=NULL;
		           rp->isUsed(false);
		           target_mrnas->Forget(rpidx);
		           target_mrnas->Delete(rpidx);
		      }
		     }
		   } //check for duplicate ref transcripts
		} //ref transcripts
		else { //-- query transfrags
		   if (m->strand=='+') { target_mrnas = &(gdata->mrnas_f); }
		     else if (m->strand=='-') { target_mrnas=&(gdata->mrnas_r); }
		        else { m->strand='.'; target_mrnas=&(gdata->umrnas); }
		   total_kept++;
		   // discard duplicate sample transfrags (but will also check for redundancy at the end)
		   if (discardDups) { //check for a redundant transfrag already loaded
			 int rpidx=-1;
			 GffObj* rp= is_TDup(m, *target_mrnas, rpidx, qDupStrict);
			 if (rp!=NULL) {
				 //always discard the shorter transfrag
				 tredundant++;
				 total_kept--;
				 if (betterTDup(rp, m)) {
					if (debug)
					   GMessage("\tQuery transcript %s discarded (duplicate of %s) (%d exons)\n",
					      m->getID(), rp->getID(), rp->exons.Count() );
					continue;
				 }
				 else {
					if (debug)
					   GMessage("\tQuery transcript %s discarded (duplicate of %s) (%d exons)\n",
					      rp->getID(), m->getID(), rp->exons.Count() );
					 ((CTData*)(rp->uptr))->mrna=NULL;
					 rp->isUsed(false);
					 target_mrnas->Forget(rpidx);
					 target_mrnas->Delete(rpidx);
				 }
			 }
		   }// redundant transfrag check
		   /* if (m->gscore==0.0)
		     m->gscore=m->exons[0]->score; //Cufflinks exon score = isoform abundance
           */

		   const char* expr = m->getAttr("FPKM");
		   if (expr!=NULL) {
		       if (expr[0]=='"') expr++;
		       fpkm=strtod(expr, NULL);
		       }
		   /* else { //backward compatibility: read RPKM if FPKM not found
		       //expr=(gtf_tracking_largeScale) ? m->getAttr("RPKM") : m->exons[0]->getAttr(m->names,"RPKM");
		       expr=m->getAttr("RPKM");
		       if (expr!=NULL) {
		           if (expr[0]=='"') expr++;
		           fpkm=strtod(expr, NULL);
		           }
		       } */
		   //const char* scov=(gtf_tracking_largeScale) ? m->getAttr("cov") : m->exons[0]->getAttr(m->names,"cov");
		   const char* scov=m->getAttr("cov");
		   if (scov!=NULL) {
		       if (scov[0]=='"') scov++;
		       cov=strtod(scov, NULL);
		       }
		   //const char* sconf_hi=(gtf_tracking_largeScale) ? m->getAttr("conf_hi") : m->exons[0]->getAttr(m->names,"conf_hi");
		   const char* stpm=m->getAttr("TPM");
		   if (stpm!=NULL){
		       if (stpm[0]=='"') stpm++;
		       tpm=strtod(stpm, NULL);
		       }
		   /*
		   //const char* sconf_lo=(gtf_tracking_largeScale) ? m->getAttr("conf_lo") : m->exons[0]->getAttr(m->names,"conf_lo");
		   const char* sconf_lo=m->getAttr("conf_lo");
		   if (sconf_lo!=NULL) {
		       if (sconf_lo[0]=='"') sconf_lo++;
		       conf_lo=strtod(sconf_lo, NULL);
		       }
		   */
		   } //query transfrags redundancy check
		target_mrnas->Add(m);
		m->isUsed(true);
		CTData* mdata=new CTData(m);
		mdata->qset=qfidx;
		gdata->tdata.Add(mdata);
		if (!is_ref_set) {
		   //mdata->dup_of=dup_by;
		//  StringTie attributes parsing
		   mdata->FPKM=fpkm;
		   mdata->cov=cov;
		   mdata->TPM=tpm;
		   //mdata->conf_lo=conf_lo;
		   }
	}//for each mrna read
	/*
	if (gtf_tracking_verbose && total_kept!=total_seen) {
		if (is_ref_set) {
          GMessage(" Kept %d ref transcripts out of %d (%d redundant discarded)\n",
        		   total_kept, total_seen, tredundant);
		}
		else {
	     GMessage(" Kept %d transfrags out of %d (%d redundant discarded)\n",
	    		 total_kept, total_seen, tredundant);
		}
	}
	*/
	return tredundant;
}

/*
bool tMatch(GffObj& a, GffObj& b, int& ovlen, bool relaxed_singleExonMatch, bool contain_only) {
	//strict intron chain match, or single-exon match
	int imax=a.exons.Count()-1;
	int jmax=b.exons.Count()-1;
	ovlen=0;
	if (imax!=jmax) return false; //different number of exons, cannot match
	if (imax==0) { //single-exon mRNAs
		if (contain_only) { //require strict boundary containment (a in b or b in a)
			//but also that at least 80% of the largest one be covered
		   if (strictMatching)
			   return (a.exons[0]->start==b.exons[0]->start &&
			   						a.exons[0]->end==b.exons[0]->end);
		   else
			   return ((a.start>=b.start && a.end<=b.end && a.covlen>=b.covlen*0.8) ||
		           (b.start>=a.start && b.end<=a.end && b.covlen>=a.covlen*0.8));
		}
		if (relaxed_singleExonMatch) { //contain_only was already tested
			return (singleExonTMatch(a,b,ovlen));
		} else {
			//same as contain_only, but stricter (at least 90% larger transcript coverage)
			if (strictMatching)
				return (a.exons[0]->start==b.exons[0]->start &&
						a.exons[0]->end==b.exons[0]->end);
			else
			   return ((a.start>=b.start && a.end<=b.end && a.covlen>=b.covlen*0.9) ||
			        (b.start>=a.start && b.end<=a.end && b.covlen>=a.covlen*0.9));
		}
	}
	if ( a.exons[imax]->start<b.exons[0]->end ||
		b.exons[jmax]->start<a.exons[0]->end )
		return false; //intron chains do not overlap at all
	//check intron overlaps
	ovlen=a.exons[0]->end-(GMAX(a.start,b.start))+1;
	ovlen+=(GMIN(a.end,b.end))-a.exons.Last()->start;
	for (int i=1;i<=imax;i++) {
		if (i<imax) ovlen+=a.exons[i]->len();
		if ((a.exons[i-1]->end!=b.exons[i-1]->end) ||
			(a.exons[i]->start!=b.exons[i]->start)) {
			return false; //intron mismatch
		}
	}
	//--- intron chain is matching ---
	if (contain_only) {//requires actual coordinate containing
		     if (strictMatching)
		    	 return (a.exons[0]->start==b.exons[0]->start &&
		    			 a.exons.Last()->end==b.exons.Last()->end);
		     else return ((a.start>=b.start && a.end<=b.end) ||
		           (b.start>=a.start && b.end<=a.end));
	}
	else return true;
}
*/

void cluster_mRNAs(GList<GffObj> & mrnas, GList<GLocus> & loci, int qfidx) {
	//mrnas sorted by start coordinate
	//and so are the loci
	//int rdisc=0;
		for (int t=0;t<mrnas.Count();t++) {
		GArray<int> mrgloci(false);
		GffObj* mrna=mrnas[t];
		int lfound=0; //count of parent loci
		/*for (int l=0;l<loci.Count();l++) {
			if (loci[l]->end<mrna->exons.First()->start) continue;
			if (loci[l]->start>mrna->exons.Last()->end) break; */
		 for (int l=loci.Count()-1;l>=0;l--) {
		   if (loci[l]->end<mrna->exons.First()->start) {
		       if (mrna->exons.First()->start-loci[l]->start > GFF_MAX_LOCUS) break;
		       continue;
		       }
		   if (loci[l]->start>mrna->exons.Last()->end) continue;
			//here we have mrna overlapping loci[l]
			if (loci[l]->add_mRNA(mrna)) {
				//a parent locus was found
				lfound++;
				mrgloci.Add(l); //locus indices added here, in decreasing order
			}
		}//loci loop
		//if (lfound<0) continue; //mrna was a ref duplicate, skip it
		if (lfound==0) {
			//create a locus with only this mRNA
 			 loci.Add(new GLocus(mrna, qfidx));
		    }
		 else if (lfound>1) {
			//more than one locus found parenting this mRNA, merge loci
		     lfound--;
			 for (int l=0;l<lfound;l++) {
				  int mlidx=mrgloci[l]; //largest indices first, so it's safe to remove
				  loci[mrgloci[lfound]]->addMerge(*loci[mlidx], mrna);
				  loci.Delete(mlidx);
			    }
		    }
	}//mrnas loop
	//if (rdisc>0) mrnas.Pack();
	//return rdisc;
}

void gatherRefLocOvls(GffObj& m, GLocus& rloc) {
	if (m.start>rloc.end || m.end<rloc.start) {
		return; //nothing to do
	}
	for (int i=0;i<rloc.mrnas.Count();i++) {
		GffObj* r=rloc.mrnas[i];
		TOvlData ovld=getOvlData(m,*r, stricterMatching);
		if (ovld.ovlcode!=0) { //has some sort of overlap with r
			((CTData*)m.uptr)->addOvl(ovld,r);
			//if (classcode_rank(olen>ovlen) { ovlen=olen; rovl=r; }
			if (ovld.ovlcode=='c' || ovld.ovlcode=='=' || ovld.ovlcode=='~') //keep match/containment for each reference transcript
				((CTData*)r->uptr)->addOvl(ovld, &m);
		}
	}//for each ref in rloc
	//GffObj** rr=&rovl;
	//char best_code=((CTData*)m.uptr)->getBestCode(rr, &ovlen);
	//return best_code;
}

int getMaxOvl(GffObj* m, GList<GffObj>& mrnas) {
	int maxovl=0;
	if (mrnas.Count()>0) {
		int qidx=qsearch_mrnas(m->end, mrnas);
		//qidx is lowest index having mrnas[qidx]->start > m->end
		// so mrnas[qidx-1]->start <= m->end
		if (qidx!=0) {
			if (qidx==-1) qidx=mrnas.Count();
			for (int i=qidx-1;i>=0;i--) {
				GffObj& t=*mrnas[i];
				if (m->start > t.end) {
					if (m->start > t.start+GFF_MAX_LOCUS)
						break; //went too far back, give up
					continue;
				}
				if (t.start>m->end) continue; //shouldn't happen
				//m overlaps t
				int ovl=m->exonOverlapLen(t);
				if (ovl>maxovl) maxovl=ovl;
			}
		}
	}
	return maxovl;
}

GffObj* gatherRefOvls(GffObj *m, GList<GLocus>& loci, int& ovlen) {
	//return the best ref overlap data for m, even when
	// no new overlaps are gathered in this call
	GffObj* r=NULL;
	if (loci.Count()>0) {
		int qidx=qsearch_loci(m->end, loci);
		//qidx is lowest index having loci[qidx]->start > m->end
		// so loci[qidx-1]->start <= m->end
		if (qidx!=0) {
			if (qidx==-1) qidx=loci.Count();
			for (int i=qidx-1;i>=0;i--) {
				GLocus& loc=*loci[i];
				if (m->start > loc.end) {
					if (m->start > loc.start+GFF_MAX_LOCUS)
						break; //went too far back, give up
					continue;
				}
				if (loc.start>m->end) continue; //shouldn't happen
				//m overlaps locus loc, check transcript overlaps
				gatherRefLocOvls(*m, loc);
			}
		}
	}
	((CTData*)m->uptr)->getBestCode(&r, &ovlen);
	return r;
}

int umrnas_assignStrand(GSeqData& seqdata, GSeqData* rdata) {
	//attempt to find the strand for seqdata.umrnas (undetermined strand q)
	//based on a) overlaps with oriented reference mRNAs if present
	//         b) overlaps with oriented mRNAs from the same input set
	// stupid complication: if there are refs overlapping this unstranded mrna on BOTH strands (?!)
	//   then we need an overlap code priority to assign the "best-overlap" strand
	int fixed=0;
	for (int j=0;j<seqdata.umrnas.Count();j++) {
		if (seqdata.umrnas[j]->strand!='.') continue;
		GffObj* m=seqdata.umrnas[j];
		if (rdata!=NULL) {
			GffObj* refovl=NULL;
			int ovlen=0;
			gatherRefOvls(m, rdata->loci_f, ovlen);
			refovl=gatherRefOvls(m, rdata->loci_r, ovlen);
			if (ovlen>0 && refovl!=NULL) {
				m->strand=refovl->strand;
				if (m->strand=='+') {
						seqdata.mrnas_f.Add(m);
						seqdata.umrnas.Forget(j);
						fixed++;
				}
					else if (m->strand=='-') {
						seqdata.mrnas_r.Add(m);
						seqdata.umrnas.Forget(j);
						fixed++;
				}
			}
		} //if rdata
	} //umrnas loop
	//refassign=fixed;
	//---- now compare to other qry transcripts that already have a strand
	int maxovl_f=0; //maximum overlap found with forward strand transcripts
	int maxovl_r=0; //maximum overlap found with reverse strand transcripts
	for (int j=0;j<seqdata.umrnas.Count();j++) {
		GffObj* m=seqdata.umrnas[j];
		if (m==NULL) continue; //already assigned
		maxovl_f=getMaxOvl(m, seqdata.mrnas_f);
		maxovl_r=getMaxOvl(m, seqdata.mrnas_r);
		if (maxovl_f>maxovl_r) {
			m->strand='+';
			seqdata.mrnas_f.Add(m);
			seqdata.umrnas.Forget(j);
			fixed++;
		} else if (maxovl_r>maxovl_f) {
			m->strand='-';
			seqdata.mrnas_r.Add(m);
			seqdata.umrnas.Forget(j);
			fixed++;
		}
	}
	if (fixed>0) seqdata.umrnas.Pack();
	//if (gtf_tracking_verbose) {
	//	GMessage(" %d out of %d (%d left, %d) unoriented transfrags were assigned a strand based on overlaps.\n", fixed, incount, seqdata.umrnas.Count(), fcount);
	//}
	return fixed;
}

//retrieve ref_data for a specific genomic sequence
GSeqData* getRefData(int gid, GList<GSeqData>& ref_data) {
	int ri=-1;
	GSeqData f(gid);
	GSeqData* r=NULL;
	if (ref_data.Found(&f,ri))
		r=ref_data[ri];
	return r;
}

void read_transcripts(FILE* f, GList<GSeqData>& seqdata,
#ifdef CUFFLINKS
  boost::crc_32_type& crc_result,
#endif
   bool keepAttrs) {
	rewind(f);
	GffReader gffr(f, true); //loading only recognizable transcript features
	gffr.showWarnings(gtf_tracking_verbose);
	//          keepAttrs    mergeCloseExons   noExonAttrs
	gffr.keepAttrs(keepAttrs, true);
	gffr.mergeCloseExons(true);
	gffr.readAll();
#ifdef CUFFLINKS
     crc_result = gffr.current_crc_result();
#endif
	//                               is_ref?    check_for_dups,
	parse_mRNAs(gffr.gflst, seqdata, false,       false);
}

int cmpGSeqByName(const pointer p1, const pointer p2) {
 return strcmp(((GSeqData*)p1)->gseq_name, ((GSeqData*)p2)->gseq_name);
}

void sort_GSeqs_byName(GList<GSeqData>& seqdata) {
  seqdata.setSorted(&cmpGSeqByName);
}

void read_mRNAs(FILE* f, GList<GSeqData>& seqdata, GList<GSeqData>* ref_data,
	         bool discardDups, int qfidx, const char* fname, bool only_multiexon) {
			 //bool intron_poking, bool keep_dups) {
	//>>>>> read all transcripts/features from a GTF/GFF3 file
	//int imrna_counter=0;
#ifdef HEAPROFILE
    if (IsHeapProfilerRunning())
      HeapProfilerDump("00");
#endif
	int loci_counter=0;
	if (ref_data==NULL) ref_data=&seqdata;
	bool isRefData=(&seqdata==ref_data);
	                          //(f, transcripts_only)
	GffReader* gffr=new GffReader(f, true); //load only transcript annotations
	gffr->showWarnings(gtf_tracking_verbose);
	//            keepAttrs=!isRefData,   mergeCloseExons   noExonAttrs=(isRefData || gtf_tracking_largeScale)
	gffr->mergeCloseExons(!noMergeCloseExons);
	const char* fext=getFileExt(fname);
	if (Gstricmp(fext, "bed")==0)
	//char* fbed=strifind(fname, ".bed");
	//if (fbed!=NULL && (size_t)(fbed-fname)>=strlen(fname)-6)
	   gffr->isBED(true);
	gffr->keepAttrs(!isRefData, isRefData || gtf_tracking_largeScale );
	gffr->readAll();
	//gffr->readAll(!isRefData,          true,        isRefData || gtf_tracking_largeScale);
	//so it will read exon attributes only for low number of Cufflinks files
#ifdef HEAPROFILE
    if (IsHeapProfilerRunning())
      HeapProfilerDump("post_readAll");
#endif
    //if (!isRefData && gtf_tracking_verbose)
	if (isRefData)
		GMessage("  %d reference transcripts loaded.\n", gffr->gflst.Count());
	else
		if (!gtf_tracking_largeScale)
			GMessage("  %d query transfrags loaded.\n", gffr->gflst.Count());
    int d=parse_mRNAs(gffr->gflst, seqdata, isRefData, discardDups, qfidx,
    		             only_multiexon);
#ifdef HEAPROFILE
    if (IsHeapProfilerRunning())
      HeapProfilerDump("post_parse_mRNAs");
#endif
	if (d>0) { //(gtf_tracking_verbose && d>0)
	  if (isRefData) GMessage("  %d duplicate reference transcripts discarded.\n",d);
	            else GMessage("  %d duplicate query transfrags discarded.\n",d);
	}
	//imrna_counter=gffr->mrnas.Count();
	delete gffr; //free the extra memory and unused GffObjs
#ifdef HEAPROFILE
    if (IsHeapProfilerRunning())
      HeapProfilerDump("post_del_gffr");
#endif

	//for each genomic sequence, cluster transcripts
	int oriented_by_overlap=0;
	int initial_unoriented=0;
	int final_unoriented=0;
	GStr bname(fname);
	GStr s;
	if (!bname.is_empty()) {
		int di=bname.rindex('.');
		if (di>0) bname.cut(di);
		int p=bname.rindex('/');
		if (p<0) p=bname.rindex('\\');
		if (p>=0) bname.remove(0,p);
	}
	FILE* fdis=NULL;
	FILE* frloci=NULL;

	for (int g=0;g<seqdata.Count();g++) {
		//find the corresponding refseqdata with the same gseq_id
		int gseq_id=seqdata[g]->get_gseqid();
		if (!isRefData) { //query data, find corresponding ref data
			GSeqData* rdata=getRefData(gseq_id, *ref_data);
			initial_unoriented+=seqdata[g]->umrnas.Count();
			if (seqdata[g]->umrnas.Count()>0) {
			    oriented_by_overlap+=umrnas_assignStrand(*seqdata[g], rdata); //, fdis);
			    final_unoriented+=seqdata[g]->umrnas.Count();
			    }
			}
		//>>>>> group mRNAs into locus-clusters (based on exon overlap)
		cluster_mRNAs(seqdata[g]->mrnas_f, seqdata[g]->loci_f, qfidx);
		cluster_mRNAs(seqdata[g]->mrnas_r, seqdata[g]->loci_r, qfidx);
		if (!isRefData) {
			cluster_mRNAs(seqdata[g]->umrnas, seqdata[g]->nloci_u, qfidx);
		}
		loci_counter+=seqdata[g]->loci_f.Count();
		loci_counter+=seqdata[g]->loci_r.Count();
//		if (refData) {
//			if (frloci==NULL) {
//				s=bname;
//				s.append(".loci.lst");
//				frloci=fopen(s.chars(), "w");
//			}
//			writeLoci(frloci, seqdata[g]->loci_f);
//			writeLoci(frloci, seqdata[g]->loci_r);
//		}//write ref loci
	}//for each genomic sequence
	if (fdis!=NULL) fclose(fdis);
	if (frloci!=NULL) fclose(frloci);
	if (initial_unoriented || final_unoriented) {
	  if (gtf_tracking_verbose) {
		        if (oriented_by_overlap>0) GMessage("  Found %d transfrags with undetermined strand (%d out of initial %d were fixed by overlaps)\n",
			     final_unoriented, oriented_by_overlap, initial_unoriented);
			else GMessage("  Found %d transfrags with undetermined strand.\n",
			     final_unoriented);
	  }
	}
	//if (fdis!=NULL) remove(s.chars()); remove 0-length file
#ifdef HEAPROFILE
    if (IsHeapProfilerRunning())
      HeapProfilerDump("post_cluster");
#endif
}

int qsearch_mrnas(uint x, GList<GffObj>& mrnas) {
	//quick search on sorted mrnas list
	//return the lowest idx where mrnas[idx]->start > x
	//---caller should make sure that mrnas.Count()>0 !
	if (mrnas[0]->start>x) return 0; //all start after x
	if (mrnas.Last()->start<x) return -1; //all start before x
	uint mstart=0;
	int mi=0;
	int idx=-1;
	int maxh=mrnas.Count()-1;
	int l=0;
	int h = maxh;
	while (l <= h) {
		mi = (l+h)>>1; //pivot index
		mstart=mrnas[mi]->start;
		if (mstart < x)  l = mi + 1; //search upper half
		else { // mstart >= x
			if (mstart == x) { //found matching start
				idx=mi;//just find the first item starting above x
				while (idx<=maxh && mrnas[idx]->start==x)
					idx++;
				return (idx>maxh) ? -1 : idx;
			}
			h = mi - 1; //search lower half
		}
	} //while there's a range of indexes to check
	idx = l; //no match found, h==l
	//--find first item starting above x
	while (idx<=maxh && mrnas[idx]->start<=x)
		idx++;
	return (idx>maxh) ? -1 : idx;
}

int qsearch_loci(uint x, GList<GLocus>& loci) {
	//quick search on sorted loci list
	//return the lowest idx where loci[idx]->start > x
	//---caller should make sure that loci.Count()>0 !
	if (loci[0]->start>x) return 0;
	if (loci.Last()->start<x) return -1;
	uint mstart=0;
	int mi=0;
	int idx=-1;
	int maxh=loci.Count()-1;
	int l=0;
	int h = maxh;
	while (l <= h) {
		mi = (l + h) >> 1;
		mstart=loci[mi]->start;
		if (mstart < x) l=mi+1;
		else {
			if (mstart == x) { //found matching coordinate here
				idx=mi;
				while (idx<=maxh && loci[idx]->start==x)
					idx++;
				return (idx>maxh) ? -1 : idx;
			}
			h=mi-1;
		}
	} //while
	idx = l;
	while (idx<=maxh && loci[idx]->start<=x) {
		idx++;
	}
	return (idx>maxh) ? -1 : idx;
}

