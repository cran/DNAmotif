#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
#include <vector>
#include <algorithm>
#include <iostream>

struct motif_struct {
    int Sl_No;
    string Motif;
};

vector<motif_struct> readMotifsFile(string inputFileName) {
    vector<motif_struct> motifs_list;
    std::ifstream inputFile(inputFileName);
    std::string line1;
    int k = 0;

    if (inputFile.is_open()) {
        while (getline(inputFile, line1)) {
            k++;
            motifs_list.push_back({k, line1});
        }
        inputFile.close();
    } else {
        Rcpp::Rcout << "Unable to open motif file" << std::endl;
    }

    return motifs_list;
}

string sequenceToString(vector<unsigned char> query, int len1) {
    unsigned char alphabet[5] = {'A', 'C', 'G', 'T', 'N'};
    std::string s;
    for (int i = 0; i < len1; ++i) {
        s.push_back(alphabet[query[i]]);
    }
    return s;
}

float dynomic_program_sim(string seq1, string seq2, int len1, int len2) {
    int len = max(len1, len2);
    vector<int> align(len, 0), tAlign(len, 0);

    for (int i = 0; i < len1; i++) {
        for (int j = 0; j < len2 - 1; j++) {
            if (seq1[i] == seq2[j])
                align[j + 1] = tAlign[j] + 1;
            else
                align[j + 1] = max(tAlign[j + 1], align[j]);
            tAlign[j + 1] = align[j + 1];
        }
    }
    return (float) align[len2 - 1] / (float) len;
}

// [[Rcpp::export(.generate_frequency_table)]]
Rcpp::DataFrame generate_frequency_table(std::string head_seq, std::string seq_id, int seq_size, int ws, float cut_off, Rcpp::StringVector motifs, int k) {
    int array_size = motifs.size();
    vector<int> arr0(array_size, 0);

    for (int i1 = 0; i1 < seq_size; i1 += ws) {
        std::string t_motif = head_seq.substr(i1, ws);
        float matched_v = 0;
        int matched_j = 0, matched_f = 0;

        for (int j = 0; j < k; j++) {
            float sim = dynomic_program_sim(Rcpp::as<std::string>(motifs[j]), t_motif, ws, ws);
            if (sim > cut_off) {
                if (sim > matched_v) {
                    matched_j = j;
                    matched_v = sim;
                    matched_f = 1;
                }
            }
        }

        if (matched_f == 1) {
            arr0[matched_j] += 1;
        } else {
            if ((int) t_motif.length() == ws) {
                k++;
                motifs.push_back(t_motif);
                arr0.push_back(1);
            }
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("motif") = motifs, Rcpp::Named("count") = arr0);
}

