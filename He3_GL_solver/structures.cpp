#include "structures.hpp"

Bound_Cond& Bound_Cond::operator= (Bound_Cond& rhs) {
    typeB = rhs.typeB;
    typeT = rhs.typeT;
    typeL = rhs.typeL;
    typeR = rhs.typeR;
    valB = rhs.valB;
    valT = rhs.valT;
    valL = rhs.valL;
    valR = rhs.valR;
    return *this;
}

ifstream& operator>> (ifstream& stream, Bound_Cond& BC) {
    stream >> BC.typeT;
    stream >> BC.valT;
    stream >> BC.typeB;
    stream >> BC.valB;
    stream >> BC.typeL;
    stream >> BC.valL;
    stream >> BC.typeR;
    stream >> BC.valR;
    return stream;
}

ostream& operator<< (ostream& out, Bound_Cond& BC) {
    out << "BC.typeT: " << BC.typeT << endl;
    out << "BC.valT:  " << BC.valT << endl;
    out << "BC.typeB: " << BC.typeB << endl;
    out << "BC.valB:  " << BC.valB << endl;
    out << "BC.typeL: " << BC.typeL << endl;
    out << "BC.valL:  " << BC.valL << endl;
    out << "BC.typeR: " << BC.typeR << endl;
    out << "BC.valR:  " << BC.valR << endl;
    return out;
}