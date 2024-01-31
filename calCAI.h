int BASE(const char x) {
    static const std::unordered_map<char, int> baseMap = {
        {'A', 1},
        {'C', 2},
        {'G', 3},
        {'U', 4}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return 0;
    }
}

std::vector<std::string> codonset={"AUU","AUA","AUC","CUA","CUC","CUG","CUU","UUA","UUG","GUU","GUA","GUC","GUG","UUU","UUC","AUG","UGU","UGC","GCA","GCC","GCG","GCU","GGU","GGC","GGA","GGG","CCU","CCC","CCA","CCG","ACU","ACC","ACA","ACG","UCU","UCC","UCA","UCG","AGU","AGC","UAU","UAC","UGG","CAA","CAG","AAU","AAC","CAU","CAC","GAA","GAG","GAU","GAC","AAA","AAG","CGU","CGC","CGA","CGG","AGA","AGG","UAA","UAG","UGA"};

std::unordered_map<int, float> codonset_CAI_;

inline int cal_CAI_index(int base_0 , int base_1, int base_2){
    return base_0*100 + base_1*10 + base_2;
}

inline void initializeCodonsetCAI(std::vector<float>& cai_vector){
    for (int i = 0; i<codonset.size(); i++){

        int tmp = cal_CAI_index(BASE(codonset[i][0]), BASE(codonset[i][1]), BASE(codonset[i][2]));

        codonset_CAI_[tmp]=log10f(cai_vector[i]);

    }    
}

float cal_CAI(std::string result_seq, std::vector<float>& cai_vector){
    initializeCodonsetCAI(cai_vector);

    float CAI=0;
    int cai_index;
    for(int i = 0; i < result_seq.size(); i++){
        cai_index = cal_CAI_index(BASE(result_seq[3*i]), BASE(result_seq[3*i+1]), BASE(result_seq[3*i+2]));
        CAI = CAI + codonset_CAI_[cai_index];
    }
    return pow(10.0, CAI/result_seq.size());
}