#include "t_classify.h"

int classcode_rank(char c) {
	switch (c) {
		case '=': return 0; //intron chain match
		case 'c': return 2; //containment, perfect partial match (transfrag < reference)
		case 'k': return 4; // reverse containment (reference < transfrag)
		case 'm': return 6; // full span overlap with all reference introns either matching or retained
		case 'n': return 8; // partial overlap transfrag with at least one intron retention
		case 'j': return 10; // multi-exon transfrag with at least one junction match
		case 'e': return 12; // single exon transfrag partially overlapping an intron of reference (possible pre-mRNA fragment)
		case 'o': return 14; // other generic exon overlap
		case 's': return 16; //"shadow" - an intron overlaps with a ref intron on the opposite strand (wrong strand mapping?)
		case 'x': return 18; // generic overlap on opposite strand (usually wrong strand mapping)
		case 'i': return 20; // intra-intron (transfrag fully contained within a reference intron)
		case 'y': return 30; // no exon overlap: ref exons fall within transfrag introns!
		case 'p': return 90; //polymerase run
		case 'r': return 92; //repeats
		case 'u': return 94; //intergenic
		case  0 : return 100;
		 default: return 96;
		}
}

bool singleExonTMatch(GffObj& m, GffObj& r, int& ovlen) {
 //if (m.exons.Count()>1 || r.exons.Count()>1..)
 GSeg mseg(m.start, m.end);
 ovlen=mseg.overlapLen(r.start,r.end);
 // fuzzy matching for single-exon transcripts:
 // overlap should be 80% of the length of the longer one
 if (m.covlen>r.covlen) {
   return ( (ovlen >= m.covlen*0.8) ||
		   (ovlen >= r.covlen*0.8 && ovlen >= m.covlen* 0.7 ));
		   //allow also some fuzzy reverse containment
 } else
   return (ovlen >= r.covlen*0.8);
}

char getOvlCode_t(GffObj& m, GffObj& r, int& ovlen) {
	ovlen=0; //total actual exonic overlap
	if (!m.overlap(r.start,r.end)) return 0;
	int jmax=r.exons.Count()-1;
	//int iovlen=0; //total m.exons overlap with ref introns
	char rcode=0;
	if (m.exons.Count()==1) { //single-exon transfrag
		GSeg mseg(m.start, m.end);
		if (jmax==0) { //also single-exon ref
			//ovlen=mseg.overlapLen(r.start,r.end);
			if (singleExonTMatch(m, r, ovlen))
				return '=';
			if (m.covlen<r.covlen)
			   { if (ovlen >= m.covlen*0.8) return 'c'; } // fuzzy containment
			else
				if (ovlen >= r.covlen*0.8 ) return 'k';   // fuzzy reverse containment
			return 'o'; //just plain overlapping
		}
		//single-exon qry overlaping multi-exon ref
		for (int j=0;j<=jmax;j++) {
			//check if it's ~contained by an exon
			int exovlen=mseg.overlapLen(r.exons[j]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (m.start>r.exons[j]->start-4 && m.end<r.exons[j]->end+4) {
					return 'c'; //close enough to be considered contained in this exon
				}
			}
			if (j==jmax) break;
			//check if it's fully contained by an intron
			if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
				return 'i';
			// check if it's a potential pre-mRNA transcript
			// (if overlaps an intron at least 10 bases)
			uint introvl=mseg.overlapLen(r.exons[j]->end+1, r.exons[j+1]->start-1);
			//iovlen+=introvl;
			if (introvl>=10 && mseg.len()>introvl+10) { rcode='e'; }
		}
		//
		if (rcode>0) return rcode;
		return 'o'; //plain overlap, uncategorized
	} //single-exon transfrag
	//-- from here on we check a multi-exon transfrag --
	int imax=m.exons.Count()-1;// imax>0 here

	if (jmax==0) { //single-exon reference overlap
		//any exon overlap?
		GSeg rseg(r.start, r.end);
		for (int i=0;i<=imax;i++) {
			//check if it's ~contained by an exon
			int exovlen=rseg.overlapLen(m.exons[i]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (r.start>m.exons[i]->start-4 && r.end<m.exons[i]->end+4) {
					return 'k'; //reference contained in this assembled exon
				}
			}
			if (i==imax) break;
			if (r.end<m.exons[i+1]->start && r.start>m.exons[i]->end)
				return 'y'; //ref contained in this transfrag intron
		}
		return 'o';
	}
	// * check if transfrag contained by a ref intron
	for (int j=0;j<jmax;j++) {
		if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
			return 'i';
	}
	//> check if m's intron chain is a subset of  r's intron chain
	if (m.exons[imax]->start<r.exons[0]->end ||
			r.exons[jmax]->start<m.exons[0]->end ) //intron chains cannot overlap
		return 'o'; //but terminal exons do!
	int i=1; //index of exon to the right of current qry intron
	int j=1; //index of exon to the right of current ref intron
	while (i<=imax && j<=jmax) { //find first intron overlap
		if (r.exons[j]->start<m.exons[i-1]->end) { j++; continue; }
		if (m.exons[i]->start<r.exons[j-1]->end) { i++; continue; }
		break; //we have an intron overlap
	}
	if (i>imax || j>jmax)
		return 'o'; //no initial intron overlap found
	//from here on we check all qry introns against ref introns
	bool junct_match=false; //true if at least a junction match is found
	//bool ichain_match=(i==1); //intron chain match, to be updated as introns are checked
	bool ichain_match=true; //if there is intron (sub-)chain match
	int imfirst=i; //index of first intron match in query
	int jmfirst=j; //index of first intron match in reference
	int imlast=-1;  //index of first intron match in query
	int jmlast=-1;  //index of first intron match in reference
	//check for intron matches
	while (i<=imax && j<=jmax) {
		uint mstart=m.exons[i-1]->end;
		uint mend=m.exons[i]->start;
		uint rstart=r.exons[j-1]->end;
		uint rend=r.exons[j]->start;
		if (rend<mstart) { j++; ichain_match=false; continue; } //no intron overlap, skipping ref intron
		if (mend<rstart) { i++; ichain_match=false; continue; } //no intron overlap, skipping qry intron
		//jovllast=j; //keep track of last overlapping intron
		//iovllast=i;
		//overlapping introns, test junction matching
		bool smatch=(mstart==rstart);
		bool ematch=(mend==rend);
		if (smatch || ematch) junct_match=true;
		if (smatch && ematch) {
			//perfect match for this intron
			if (ichain_match) { imlast=i; jmlast=j; }
			i++; j++;
			continue;
		}//intron match
		//-- not an intron match:
		ichain_match=false;
		if (mend>rend) j++; else i++;
	} //while checking intron overlaps

	if (ichain_match) { //intron sub-chain match
		if (imfirst==1 && imlast==imax) { // qry full intron chain match
			if (jmfirst==1 && jmlast==jmax) return '='; //identical intron chains
			// -- qry intron chain is shorter than ref intron chain --
			int l_iovh=0;   // overhang of leftmost q exon left boundary beyond the end of ref intron to the left
			int r_iovh=0;   // same type of overhang through the ref intron on the right
			if (jmfirst>1 && r.exons[jmfirst-1]->start>m.start)
				l_iovh = r.exons[jmfirst-1]->start - m.start;
			if (jmlast<jmax && m.end > r.exons[jmlast]->end)
				r_iovh = m.end - r.exons[jmlast]->end;
			if (l_iovh<4 && r_iovh<4) return 'c';
		} else if ((jmfirst==1 && jmlast==jmax)) {//ref full intron chain match
			//check if the reference i-chain is contained in qry i-chain
			int l_jovh=0;   // overhang of leftmost q exon left boundary beyond the end of ref intron to the left
			int r_jovh=0;   // same type of overhang through the ref intron on the right
			if (imfirst>1 && m.exons[imfirst-1]->start>r.start)
				l_jovh = m.exons[imfirst-1]->start - r.start;
			if (imlast<imax && r.end > m.exons[imlast]->end)
				r_jovh = r.end - m.exons[imlast]->end;
			if (l_jovh<4 && r_jovh<4) return 'k'; //reverse containment
		}
		return 'j';
	}
	if (junct_match) return 'j';
	//we could have 'o' or 'y' here
	//any real exon overlaps?
	ovlen=m.exonOverlapLen(r);
	if (ovlen>4) return 'o';
	return 'y'; //all reference exons are within transfrag introns!
}

//from gffcompare.cpp
char getOvlCode(GffObj& m, GffObj& r, int& ovlen) {
	ovlen=0; //total actual exonic overlap
	if (!m.overlap(r.start,r.end)) return 0;
	int jmax=r.exons.Count()-1;
	//int iovlen=0; //total m.exons overlap with ref introns
	char rcode=0;
	if (m.exons.Count()==1) { //single-exon transfrag
		GSeg mseg(m.start, m.end);
		if (jmax==0) { //also single-exon ref
			//ovlen=mseg.overlapLen(r.start,r.end);
			if (singleExonTMatch(m, r, ovlen))
				return '=';
			if (m.covlen<r.covlen)
			   { if (ovlen >= m.covlen*0.8) return 'c'; } // fuzzy containment
			else
				if (ovlen >= r.covlen*0.8 ) return 'k';   // fuzzy reverse containment
			return 'o'; //just plain overlapping
		}
		//single-exon qry overlaping multi-exon ref
		//TODO: special case: single-exon transfrag fully covering an intron (code 'f')
		for (int j=0;j<=jmax;j++) {
			//check if it's ~contained by an exon
			int exovlen=mseg.overlapLen(r.exons[j]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (m.start>r.exons[j]->start-4 && m.end<r.exons[j]->end+4) {
					return 'c'; //close enough to be considered contained in this exon
				}
			}
			if (j==jmax) break;
			//check if it's fully contained by an intron
			if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
				return 'i';
			// check if it's a potential pre-mRNA transcript
			// (if overlaps this intron at least 10 bases)
			uint introvl=mseg.overlapLen(r.exons[j]->end+1, r.exons[j+1]->start-1);
			//iovlen+=introvl;
			if (introvl>=10 && mseg.len()>introvl+10) { rcode='e'; }
		}
		//
		if (rcode>0) return rcode;
		return 'o'; //plain overlap, uncategorized
	} //single-exon transfrag
	//-- from here on we check a multi-exon transfrag --
	int imax=m.exons.Count()-1;// imax>0 here

	if (jmax==0) { //single-exon reference overlap
		//any exon overlap?
		GSeg rseg(r.start, r.end);
		for (int i=0;i<=imax;i++) {
			//check if it's ~contained by an exon
			int exovlen=rseg.overlapLen(m.exons[i]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (r.start>m.exons[i]->start-4 && r.end<m.exons[i]->end+4) {
					return 'k'; //reference contained in this assembled exon
				}
			}
			if (i==imax) break;
			if (r.end<m.exons[i+1]->start && r.start>m.exons[i]->end)
				return 'y'; //ref contained in this transfrag intron
		}
		return 'o';
	}
	// * check if transfrag contained by a ref intron
	for (int j=0;j<jmax;j++) {
		if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
			return 'i';
	}
	//> check if m's intron chain is a subset of  r's intron chain
	if (m.exons[imax]->start<r.exons[0]->end ||
			r.exons[jmax]->start<m.exons[0]->end ) //intron chains cannot overlap
		return 'o'; //but terminal exons do!
	int i=1; //index of exon to the right of current qry intron
	int j=1; //index of exon to the right of current ref intron
	while (i<=imax && j<=jmax) { //find first intron overlap
		if (r.exons[j]->start<m.exons[i-1]->end) { j++; continue; }
		if (m.exons[i]->start<r.exons[j-1]->end) { i++; continue; }
		break; //we have an intron overlap
	}
	if (i>imax || j>jmax)
		return 'o'; //no initial intron overlap found
	//from here on we check all qry introns against ref introns
	bool junct_match=false; //true if at least a junction match is found
	//bool ichain_match=(i==1); //intron chain match, to be updated as introns are checked
	bool ichain_match=true; //if there is intron (sub-)chain match
	int imfirst=i; //index of first intron match in query
	int jmfirst=j; //index of first intron match in reference
	int imlast=-1;  //index of first intron match in query
	int jmlast=-1;  //index of first intron match in reference
	//check for intron matches
	//TODO: keep track of intron retention! (code 'n')
	while (i<=imax && j<=jmax) {
		uint mstart=m.exons[i-1]->end;
		uint mend=m.exons[i]->start;
		uint rstart=r.exons[j-1]->end;
		uint rend=r.exons[j]->start;
		if (rend<mstart) { j++; ichain_match=false; continue; } //no intron overlap, skipping ref intron
		if (mend<rstart) { i++; ichain_match=false; continue; } //no intron overlap, skipping qry intron
		//jovllast=j; //keep track of last overlapping intron
		//iovllast=i;
		//overlapping introns, test junction matching
		bool smatch=(mstart==rstart);
		bool ematch=(mend==rend);
		if (smatch || ematch) junct_match=true;
		if (smatch && ematch) {
			//perfect match for this intron
			if (ichain_match) { imlast=i; jmlast=j; }
			i++; j++;
			continue;
		}//intron match
		//-- not an intron match:
		ichain_match=false;
		if (mend>rend) j++; else i++;
	} //while checking intron overlaps

	if (ichain_match) { //intron sub-chain match
		if (imfirst==1 && imlast==imax) { // qry full intron chain match
			if (jmfirst==1 && jmlast==jmax) return '='; //identical intron chains
			// -- qry intron chain is shorter than ref intron chain --
			int l_iovh=0;   // overhang of leftmost q exon left boundary beyond the end of ref intron to the left
			int r_iovh=0;   // same type of overhang through the ref intron on the right
			if (jmfirst>1 && r.exons[jmfirst-1]->start>m.start)
				l_iovh = r.exons[jmfirst-1]->start - m.start;
			if (jmlast<jmax && m.end > r.exons[jmlast]->end)
				r_iovh = m.end - r.exons[jmlast]->end;
			if (l_iovh<4 && r_iovh<4) return 'c';
		} else if ((jmfirst==1 && jmlast==jmax)) {//ref full intron chain match
			//check if the reference i-chain is contained in qry i-chain
			int l_jovh=0;   // overhang of leftmost q exon left boundary beyond the end of ref intron to the left
			int r_jovh=0;   // same type of overhang through the ref intron on the right
			if (imfirst>1 && m.exons[imfirst-1]->start>r.start)
				l_jovh = m.exons[imfirst-1]->start - r.start;
			if (imlast<imax && r.end > m.exons[imlast]->end)
				r_jovh = r.end - m.exons[imlast]->end;
			if (l_jovh<4 && r_jovh<4) return 'k'; //reverse containment
		}
		return 'j';
	}
	if (junct_match) return 'j';
	//we could have 'o' or 'y' here
	//any real exon overlaps?
	ovlen=m.exonOverlapLen(r);
	if (ovlen>4) return 'o';
	return 'y'; //all reference exons are within transfrag introns!
}
