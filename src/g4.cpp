#include "g4.h"

void makejson(string name, nt value, stringstream &out) {
    out << "\"" << name << "\": " << value;
}
void makejson(string name, short value, stringstream &out) {
    out << "\"" << name << "\": " << value;
}
void makejson(string name, string value, stringstream &out) {
    out << "\"" << name << "\": \"" << value << "\"";
}

string G4 :: toString(bool print_overlaps=true) {
    stringstream out;
    out << "{";
    makejson("start", start, out); out << ",";
    makejson("tetrad1", tetrad1, out); out << ",";
    makejson("tetrad2", tetrad2, out); out << ",";
    makejson("tetrad3", tetrad3, out); out << ",";
    makejson("tetrad4", tetrad4, out); out << ",";
    makejson("y1", y1, out); out << ",";
    makejson("y2", y2, out); out << ",";
    makejson("y3", y3, out); out << ",";
    makejson("tetrads", tetrads, out); out << ",";
    makejson("length", length, out); out << ",";
    makejson("gscore", gscore, out); out << ",";
    makejson("sequence", sequence, out);
    if ( print_overlaps ) {
        out << ",";
        out << "\"overlaps\":  [";
        int i = 0;
        for ( vector<G4>::iterator git = overlaps.begin(); git != overlaps.end(); ++git) {
            out << git->toString(false);
            if ( ++i != overlaps.size()) {
                out << ",";
            }
        }
        out << "]";
    }
    out << "}";
    return out.str();
}
