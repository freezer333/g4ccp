
#include <iostream>
#include <string>
#include <queue>
#include <cmath>
#include <sstream>
#include <ctime>
using namespace std;


typedef unsigned short nt;
class G4Candidate;
class G4;



inline int maximumLength (int numTetrads) {
  // Note the QGRS Mapper requires the maximum length to be set
  // externally.  The problem with this is that the 30 was originally
  // meant to allow GG, and longer motifs should use 45.  This algorithm
  // uses dynamic assignment of maximum length.
  return (numTetrads < 3)  ? 30 : 45;
}

class G4Candidate {
public:
    G4Candidate(string sequence, short tetrads, nt start_pos) {
        this->y1 = -1;
        this->y2 = -1;
        this->y3 = -1;
        this->sequence = sequence;
        this->numTetrads = tetrads;
        this->start = start_pos;
        this->tstring = "";
        for ( int i = 0; i < tetrads; i++ ) tstring.append("G");
        this->maxLength = maximumLength(tetrads);
        //cout << "Created candidate with " << numTetrads << " -> " << tstring << endl;
    }
    G4Candidate (short tetrads, short y1, short y2, short y3) {
      this->y1 = y1;
      this->y2 = y2;
      this->y3 = y3;
      this->numTetrads = tetrads;
      this->maxLength = maximumLength(this->numTetrads);
      this->tstring = "";
      for ( int i = 0; i < tetrads; i++ ) tstring.append("G");
      string generated_loop1 = "";
      string generated_loop2 = "";
      string generated_loop3 = "";
      for ( int i = 0; i < y1; i++ ) generated_loop1.append("C");
      for ( int i = 0; i < y2; i++ ) generated_loop2.append("C");
      for ( int i = 0; i < y3; i++ ) generated_loop3.append("C");
      this->sequence = tstring + generated_loop1 + tstring + generated_loop2 + tstring + generated_loop3 + tstring;
    }
    short y1 ;
    short y2 ;
    short y3 ;
    string sequence;
    short numTetrads;
    nt start ;
    string tstring;
    short maxLength;
    string toString() {
        stringstream sstr;
        sstr << "Start = " << start << ", numtetrads = " << numTetrads;
        return sstr.str();
    }
    short score() {
        double gavg = (abs(y1-y2) + abs(y2-y3) + abs(y1-y3))/3.0;
        return floor(gmax() - gavg + gmax() * (numTetrads-2));
    }
    short score_mapper(short max) {
        double gavg = (abs(y1-y2) + abs(y2-y3) + abs(y1-y3))/2.0;
        short gmax = max - 9;
        return floor(gmax - gavg + gmax * (numTetrads-2));
    }
    short gmax(){
      // Unlike the QGRS Mapper, this value is dynamically determined based 
      // on the actual tetrad length of this motif.  This is logical, however
      // it is inconsistent with what the mapper will output.
      return maxLength - (numTetrads * 4 + 1);
    }

    short length() {
       return 4 * numTetrads + y1 + y2 + y3;
    }

    nt t1(){
        return start;
    }

    nt t2() {
        return t1() + numTetrads + y1;
    }

    nt t3() {
        return t2() + numTetrads + y2;
    }

    nt t4() {
        return t3() + numTetrads + y3;
    }

    nt cursor() {
        if (y1 < 0 ) return t1() + numTetrads;
        else if (y2 < 0) return t2() + numTetrads;
        else if (y3 < 0) return t3() + numTetrads;
        else return -1;
    }

    short partialLength() {
        short length = numTetrads * 4;
        // add the minimum loops left
        if (y1 >= 0 && y2 <0 ) {
            // only first loop is known
            if (y1 == 0)
                // other two must be at least 2
                length += 2;
            else
                length += 1;
        }
        else if (y2 >= 0 && y3 <0) {
            //first two loop lengths are known
            if (y1 == 0 || y2 == 0 ) {
                length+= 1;
            }
        }
        // add the current loops
        if (y1 > 0 ) length += y1;
        if (y2 > 0 ) length += y2;
        if (y3 > 0 ) length += y3;
        return length;
    }

    short minAcceptableLoopLength(){
        if (y1 == 0 || y2 == 0 || y3 == 0) return 1;
        else return 0;
    }

    bool complete() {
        if (y1 < 0 || y2 < 0 || y3 < 0 ) return false;
        return true;
    }

    bool viable(int min_score) {
        if (score() < min_score )
            return false;
        if (length() > maxLength )
            return false;

        // only one loop is allowed to have a 0 length
        short count = 0;
        if (y1 < 1) count+= 1;
        if (y2 < 1) count+= 1;
        if (y3 < 1) count+= 1;
        return count < 2;
    }

    void findLoopLengthsFrom(queue<int> & ys, int i) {
        int p = i;
        bool done = false;
        while (!done) {
            p = sequence.find(tstring, p);
            if (p < (start+maxLength+1) && p >= 0) {
                int y = p - i;
                if (y >= minAcceptableLoopLength() && (p-start+tstring.length()-1) < maxLength) {
                    ys.push(y);
                }
                else done = true;
            }
            else done = true;
            p += 1;
        }
    }


    void expand(queue<G4Candidate> &cands) {
        queue<int> ys;
        findLoopLengthsFrom(ys, cursor());
        while (!ys.empty() ){
            int y = ys.front();
            ys.pop();
            G4Candidate cand(sequence, numTetrads, start);
            cand.y1 = y1;
            cand.y2 = y2;
            cand.y3 = y3;
            if (y1 < 0 ) cand.y1 = y;
            else if ( y2 < 0 ) cand.y2 = y;
            else if ( y3 < 0 ) cand.y3 = y;

            if (cand.partialLength() <= cand.maxLength )
                cands.push(cand);
        }
    }
};


class G4 {
public:
    G4() {}
    G4(G4Candidate &candidate) {
        start = candidate.start;
        tetrads = candidate.numTetrads;
        tetrad1 = candidate.t1();
        tetrad2 = candidate.t2();
        tetrad3 = candidate.t3();
        tetrad4 = candidate.t4();
        y1 = candidate.y1;
        y2 = candidate.y2;
        y3 = candidate.y3;
        length = candidate.length();
        gscore = candidate.score();
        sequence = candidate.sequence.substr(candidate.start, candidate.length());
    }

    bool isequal(const G4 & other) {
        if ( start != other.start ) return false;
        if ( tetrads != other.tetrads ) return false;
        if ( y1 != other.y1) return false;
        if ( y2 != other.y2) return false;
        if ( y3 != other.y3) return false;
        return true;
    }

    string toString(bool);

    nt start;
    nt tetrad1;
    nt tetrad2;
    nt tetrad3;
    nt tetrad4;
    short y1;
    short y2;
    short y3;
    short tetrads;
    short length;
    short gscore;
    string sequence;
    vector<G4> overlaps;
};
