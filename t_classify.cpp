#include "t_classify.h"

int classcode_rank(char c) {
	switch (c) {
		case '=': return 0; //intron chain match
		case 'c': return 2; //containment, perfect partial match (transfrag < reference)
		case 'k': return 6; // reverse containment (reference < transfrag)
		case 'm': return 6; // full span overlap with all reference introns either matching or retained
		case 'n': return 6; // partial overlap transfrag with at least one intron retention
		case 'j': return 6; // multi-exon transfrag with at least one junction match
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

//formerly in gffcompare
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
		//-- single-exon qry overlaping multi-exon ref
		//check full pre-mRNA case (all introns retained): code 'm'

		if (m.start<=r.exons[0]->end && m.end>=r.exons[jmax]->start)
			return 'm';

		for (int j=0;j<=jmax;j++) {
			//check if it's ~contained by an exon
			int exovlen=mseg.overlapLen(r.exons[j]);
			if (exovlen>0) {
				ovlen+=exovlen;
				if (m.start>r.exons[j]->start-4 && m.end<r.exons[j]->end+4) {
					return 'c'; //close enough to be considered contained in this exon
				}
			}
			if (j==jmax) break; //last exon here, no intron to check
			//check if it fully covers an intron (retained intron)
			if (m.start<r.exons[j]->end && m.end>r.exons[j+1]->start)
				return 'n';
			//check if it's fully contained by an intron
			if (m.end<r.exons[j+1]->start && m.start>r.exons[j]->end)
				return 'i';
			// check if it's a potential pre-mRNA transcript
			// (if overlaps this intron at least 10 bases)
			uint introvl=mseg.overlapLen(r.exons[j]->end+1, r.exons[j+1]->start-1);
			//iovlen+=introvl;
			if (introvl>=10 && mseg.len()>introvl+10) { rcode='e'; }
		} //for each ref exon
		if (rcode>0) return rcode;
		return 'o'; //plain overlap, uncategorized
	} //single-exon transfrag
	//-- multi-exon transfrag --
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
	if (m.exons[imax]->start<r.exons[0]->end) {
		//qry intron chain ends before ref intron chain starts
		//check if last qry exon plugs the 1st ref intron
		if (m.exons[imax]->start<=r.exons[0]->end &&
			m.exons[imax]->end>=r.exons[1]->start) return 'n';
		return 'o'; //only terminal exons overlap
	}
	else if (r.exons[jmax]->start<m.exons[0]->end) {
		//qry intron chain starts after ref intron chain ends
		//check if first qry exon plugs the last ref intron
		if (m.exons[0]->start<=r.exons[jmax-1]->end &&
			m.exons[0]->end>=r.exons[jmax]->start) return 'n';
		return 'o'; //only terminal exons overlap
	}
	//check intron chain overlap (match, containment, intron retention etc.)
	int i=1; //index of exon to the right of current qry intron
	int j=1; //index of exon to the right of current ref intron
	bool intron_conflict=false; //used for checking for retained introns
	//from here on we check all qry introns against ref introns
	bool junct_match=false; //true if at least a junction match is found
	bool ichain_match=true; //if there is intron (sub-)chain match, to be updated by any mismatch
	bool intron_ovl=false; //if any intron overlap is found
	bool intron_retention=false; //if any ref intron is covered by a qry exon
	int imfirst=0; //index of first intron match in query (valid>0)
	int jmfirst=0; //index of first intron match in reference (valid>0)
	int imlast=0;  //index of first intron match in query
	int jmlast=0;  //index of first intron match in reference
	//check for intron matches
	while (i<=imax && j<=jmax) {
		uint mstart=m.exons[i-1]->end;
		uint mend=m.exons[i]->start;
		uint rstart=r.exons[j-1]->end;
		uint rend=r.exons[j]->start;
		if (rend<mstart) { //qry intron starts after ref intron ends
			if (!intron_conflict && r.exons[j]->overlap(mstart+1, mend-1))
				intron_conflict=true;
			if (!intron_retention && rstart>=m.exons[i-1]->start)
				intron_retention=true;
			if (intron_ovl) ichain_match=false;
			j++;
			continue;
		} //no intron overlap, skipping ref intron
		if (rstart>mend) { //qry intron ends before ref intron starts
			//if qry intron overlaps the exon on the left, we have an intron conflict
			if (!intron_conflict && r.exons[j-1]->overlap(mstart+1, mend-1))
				intron_conflict=true;
			if (!intron_retention && rend<=m.exons[i]->end)
				intron_retention=true;
			if (intron_ovl) ichain_match=false;
			i++;
			continue;
		} //no intron overlap, skipping qry intron
		intron_ovl=true;
		//overlapping introns, test junction matching
		bool smatch=(mstart==rstart); //TODO: what if the introns differ just by 2 bases at one end?
		bool ematch=(mend==rend);
		if (smatch || ematch) junct_match=true;
		if (smatch && ematch) {
			//perfect match for this intron
			if (ichain_match) { //chain matching still possible
			  if (jmfirst==0) jmfirst=j;
			  if (imfirst==0) imfirst=i;
			  imlast=i;
			  jmlast=j;
			}
			i++; j++;
			continue;
		}
		//intron overlapping but with at least a junction mismatch
		intron_conflict=true;
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
	}
	//'=', 'c' and 'k' where checked and assigned, check for 'm' and 'n' before falling back to 'j'
	if (!intron_conflict && (m.start<=r.exons[0]->end && m.end>=r.exons[jmax]->start)) {
			return 'm';
	}
	if (intron_retention) return 'n';
	if (junct_match) return 'j';
	//we could have 'o' or 'y' here
	//any real exon overlaps?
	ovlen=m.exonOverlapLen(r);
	if (ovlen>4) return 'o';
	return 'y'; //all reference exons are within transfrag introns!
}
