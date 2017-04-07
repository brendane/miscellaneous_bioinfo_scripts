// compile with std=c++11
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<utility>
#include<vector>

#include<math.h>
#include<stdlib.h>

using std::cerr;
using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::pair;
using std::sort;
using std::string;
using std::stringstream;
using std::vector;

struct Record {
    float maf;
    unsigned pos;
    string chr;
    vector<float> data;
    string rs;
    unsigned nChrs;
    long assignment; 

    Record():assignment(-1){}
};

// To sort by MAF in descending order
// UPDATE: now also sorts by number of genotyped strains in
// descending order after sorting by MAF.
//
bool recordRevSort(const Record &a, const Record &b) {
    return (a.maf > b.maf) && (a.nChrs > b.nChrs);
}
bool recordRevGenoSort(const Record &a, const Record &b) {
    return (a.nChrs > b.nChrs);
}

void getFreqs(const Record &x, const Record &y,
        float &pa, float &pb, float &pab) {
    unsigned ab = 0;
    unsigned a = 0;
    unsigned b = 0;
    unsigned c = 0;
    for(unsigned i=0; i < x.data.size(); i++) {
        if(std::isnan(x.data[i]) || std::isnan(y.data[i])) {
            continue;
        }
        c++;
        if(x.data[i] > 0) a++;
        if(y.data[i] > 0) b++;
        if(x.data[i] && y.data[i] > 0) ab++;
    }
    pa  = ((float)a)/c;
    pb  = ((float)b)/c;
    pab = ((float)ab)/c;
}

// Returns R2 or -1 if it is not possible for R2 to be >= thresh.
// Formula for determining possible R2 comes from Wray 2005.
float calcR2(const Record &x, const Record &y, float thresh=0) {
    float pa, pb, pab;
    if(y.maf < ( (thresh * x.maf) / (1 - x.maf * (1 - thresh)) ) ||
            y.maf >  ( x.maf / (x.maf * (1 - thresh) + thresh) )) {
        return -1.0;
    } else {
        getFreqs(x, y, pa, pb, pab);
        return pow((pab - (pa * pb)), 2) / (pa * (1-pa) * pb * (1-pb));
    }
}

int main(int argc, char* argv[]) {
    if(argc < 4) {
        cerr << "Error: wrong number of arguments" << endl;
        cerr << "r2_groups <r2 threshold> <missing threshold> <input file>" << endl;
        return 1;
    }

    bool use_rs = false;
    if(argc == 5) {
        use_rs = (bool) atoi(argv[3]);
    }

    float threshold = atof(argv[1]);
    int miss_thresh = atoi(argv[2]);

    vector<Record> data;
    ifstream instream(argv[argc-1]);
    string line;
    float aa;
    bool multiAllele;
    unsigned nChrs;

    unsigned i = 0;
    while(std::getline(instream, line)) {
        multiAllele = false;
        nChrs = 0;
        if(i != 0) {
            Record rec;
            rec.rs = "";
            unsigned c = 0;
            float refCount = 0.;
            string column;
            stringstream lineStream(line);
            bool na = false;
            while(getline(lineStream, column, '\t')) {
                na = false;
                if(c == 0) {
                    rec.chr = column;
                } else if(c == 1) {
                    rec.pos = atoi(column.c_str());
                } else if(c == 2 && use_rs) {
                    rec.rs = column;
                } else {
                    if(column == "NA" || column == "NaN" || column == "nan") {
                        na = true;
                        aa = std::nanf("");
                    } else {
                        aa = floor(atof(column.c_str()));
                    }
                    if(!na && aa > 1) {
                        multiAllele = true;
                        break;
                    }
                    if(!na) {
                        nChrs++;
                        refCount += aa;
                    }
                    rec.data.push_back(aa);
                }
                c++;
            }
            if(multiAllele) {
                continue;
            } else {
                // No variation at this site
                if(refCount == 0.0 || refCount == nChrs)
                    continue;
                rec.maf = refCount / nChrs;
                if(rec.maf > 0.5)
                    rec.maf = 1.0 - rec.maf;
                rec.nChrs = nChrs;
                data.push_back(rec);
            }
        }
        i++;
    }

    // Sort data by maf and genotyping rate
    sort(data.begin(), data.end(), recordRevSort);

    // Calculate R2 and group the variants
    long assign = 0;
    float r2;
    for(unsigned j = 0; j < data.size()-1; j++) {
        if(data[j].assignment != -1) continue;
        assign++;
        data[j].assignment = assign;
        for(unsigned k = j+1; k < data.size(); k++) {
            if(data[k].assignment != -1) continue;
            r2 = calcR2(data[j], data[k], threshold);
            if(r2 == -1.0) break;
            if(r2 >= threshold)
                data[k].assignment = assign;
        }
    }
    if(data[data.size()-1].assignment == -1)
        data[data.size()-1].assignment = assign + 1;

    // Write the output
    if(use_rs) {
        cout << "chrom\tpos\trs\tmaf\tn\tgroup" << endl;
        for(unsigned j = 0; j < data.size(); j++) {
            cout << data[j].chr << "\t" << data[j].pos << "\t" <<
                data[j].rs << "\t" <<
                data[j].maf << "\t" << data[j].nChrs << "\t" <<
                data[j].assignment << endl;
        }
    } else {
        cout << "chrom\tpos\tmaf\tn\tgroup" << endl;
        for(unsigned j = 0; j < data.size(); j++) {
            cout << data[j].chr << "\t" << data[j].pos << "\t" <<
                data[j].maf << "\t" << data[j].nChrs << "\t" <<
                data[j].assignment << endl;
        }
    }
}
