#define VALUE_MIN std::numeric_limits<int>::lowest()

#define SINGLE_MAX_LEN 30
#define MIN_STACKPAIR_GAP 5
#define MIN_MULTI_GAP 10
#define NOTON 16 // NUM_OF_TYPE_OF_NUCS
#define HAIRPIN_GAP 3 // NUM_OF_SMALLEST_DIS_CAN_PAIR

//Amino acids label
#define normal_ami -1
#define Arg 0
#define Leu 1
#define Ser 2

#define base_number 16
#define square_base_number 256


//global variables
float LDCONST = 100;
float lambda;
int beamsize;
bool _allowed_pairs[NOTON][NOTON];
bool xcon_table[3][2][NOTON][NOTON];
std::unordered_map<char, std::unordered_map<int, int>> CodonSetCAIMap;

// declaration of functions
int BASE(const std::string& x);

std::unordered_map<std::string, std::string>  Amino_label = {
        {"A","GCN"},{"C","UGY"},{"D","GAY"},{"E","GAR"},{"F","UUY"},
        {"G","GGN"},{"H","CAY"},{"I","AUH"},{"K","AAR"},{"L","YVn"},
        {"M","AUG"},{"N","AAY"},{"P","CCN"},{"Q","CAR"},{"R","MOn"},
        {"S","WSz"},{"T","ACN"},{"V","GUN"},{"W","UGG"},{"Y","UAY"},
        {"*","NNN"}
    };

// vector<std::vector<string>> CodonSet={{"I","AUU"},{"I","AUA"},{"I","AUC"},{"L","CUA"},{"L","CUC"},{"L","CUG"},{"L","CUU"},{"L","UUA"},{"L","UUG"},{"V","GUU"},"GUA","GUC","GUG","UUU","UUC","AUG","UGU","UGC","GCA","GCC","GCG","GCU","GGU","GGC","GGA","GGG","CCU","CCC","CCA","CCG","ACU","ACC","ACA","ACG","UCU","UCC","UCA","UCG","AGU","AGC","UAU","UAC","UGG","CAA","CAG","AAU","AAC","CAU","CAC","GAA","GAG","GAU","GAC","AAA","AAG","CGU","CGC","CGA","CGG","AGA","AGG","UAA","UAG","UGA"};
std::vector<std::vector<std::string>> CodonSet=
{
    {"I","AUU"},{"I","AUA"},{"I","AUC"},
    {"L","CUA"},{"L","CUC"},{"L","CUG"},{"L","CUU"},{"L","UUA"},{"L","UUG"},
    {"V","GUU"},{"V","GUA"},{"V","GUC"},{"V","GUG"},
    {"F","UUU"},{"F","UUC"},
    {"M","AUG"},
    {"C","UGU"},{"C","UGC"},
    {"A","GCA"},{"A","GCC"},{"A","GCG"},{"A","GCU"},
    {"G","GGU"},{"G","GGC"},{"G","GGA"},{"G","GGG"},
    {"P","CCU"},{"P","CCC"},{"P","CCA"},{"P","CCG"},
    {"T","ACU"},{"T","ACC"},{"T","ACA"},{"T","ACG"},
    {"S","UCU"},{"S","UCC"},{"S","UCA"},{"S","UCG"},{"S","AGU"},{"S","AGC"},
    {"Y","UAU"},{"Y","UAC"},
    {"W","UGG"},
    {"Q","CAA"},{"Q","CAG"},
    {"N","AAU"},{"N","AAC"},
    {"H","CAU"},{"H","CAC"},
    {"E","GAA"},{"E","GAG"},
    {"D","GAU"},{"D","GAC"},
    {"K","AAA"},{"K","AAG"},
    {"R","CGU"},{"R","CGC"},{"R","CGA"},{"R","CGG"},{"R","AGA"},{"R","AGG"},
    {"Stop","UAA"},{"Stop","UAG"},{"Stop","UGA"}
};

std::unordered_map<std::string, std::string>  re_Amino_label = 
{
    {"GCU","A"}, {"GCC","A"}, {"GCA","A"}, {"GCG","A"},
    {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
    {"AAU","N"}, {"AAC","N"},
    {"GAU","D"}, {"GAC","D"},
    {"UGU","C"}, {"UGC","C"},
    {"CAA","Q"}, {"CAG","Q"},
    {"GAA","E"}, {"GAG","E"},
    {"GGU","G"}, {"GGC","G"}, {"GGA","G"}, {"GGG","G"},
    {"CAU","H"}, {"CAC","H"},
    {"AUU","I"}, {"AUC","I"}, {"AUA","I"},
    {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
    {"AAA","K"}, {"AAG","K"},
    {"AUG","M"},
    {"UUU","F"}, {"UUC","F"},
    {"CCU","P"}, {"CCC","P"}, {"CCA","P"}, {"CCG","P"},
    {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
    {"ACU","T"}, {"ACC","T"}, {"ACA","T"}, {"ACG","T"},
    {"UGG","W"},
    {"UAU","Y"}, {"UAC","Y"},
    {"GUU","V"}, {"GUC","V"}, {"GUA","V"}, {"GUG","V"}

};

std::unordered_map<char, std::vector<int>> Base_table={
    {'N',{BASE("A"),BASE("C"),BASE("G"),BASE("U")}},
    {'z',{BASE("A"),BASE("C"),BASE("G"),BASE("U"),BASE("U1"),BASE("C1")}},
    {'n',{BASE("A1"),BASE("G1"),BASE("A2"),BASE("G2"),BASE("C"),BASE("U")}},
    {'A',{BASE("A")}},
    {'U',{BASE("U")}},
    {'G',{BASE("G")}},
    {'C',{BASE("C")}},
    {'W',{BASE("A"),BASE("U")}},
    {'S',{BASE("G"),BASE("C")}},
    {'M',{BASE("A"),BASE("C")}},
    {'R',{BASE("A"),BASE("G")}},
    {'Y',{BASE("C"),BASE("U")}},
    {'H',{BASE("A"),BASE("C"),BASE("U")}},
    {'V',{BASE("U1"),BASE("U_AG2"),BASE("U_AG1")}},
    {'O',{BASE("G_CU"),BASE("G_AG2"),BASE("G_AG1")}}
};

class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }

        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

enum Manner {
    NONE=0,
    MANNER_NONEtoN,
    MANNER_NONEtoF,
    MANNER_NONEtoS,
    MANNER_N_EtoN,
    MANNER_NtoC,
    MANNER_S_EtoS,
    MANNER_CtoM1,
    MANNER_M1_CtoM2,
    MANNER_M1_EtoM1,
    MANNER_M2toM1,
    MANNER_M2toMulti,
    MANNER_S_M2toMulti,
    MANNER_Multi_EtoMulti,
    MANNER_MultitoC,
    MANNER_F_EtoF,
    MANNER_CtoF,
    MANNER_F_CtoF,
    MANNER_C_StoCS,
    MANNER_S_CStoC,
    MANNER_CStoC,
    MANNER_CtoC,
    MANNER_S_CtoC
};

std::string mannerToString(Manner manner) {
    switch (manner) {
        case NONE:
            return "NONE";
        case MANNER_NONEtoN:
            return "MANNER_NONEtoN";
        case MANNER_NONEtoF:
            return "MANNER_NONEtoF";
        case MANNER_NONEtoS:
            return "MANNER_NONEtoS";
        case MANNER_N_EtoN:
            return "MANNER_N_EtoN";
        case MANNER_NtoC:
            return "MANNER_NtoC";
        case MANNER_S_EtoS:
            return "MANNER_S_EtoS";
        case MANNER_CtoC:
            return "MANNER_CtoC";
        case MANNER_S_CtoC:
            return "MANNER_S_CtoC";
        case MANNER_CtoM1:
            return "MANNER_CtoM1";
        case MANNER_M1_CtoM2:
            return "MANNER_M1_CtoM2";
        case MANNER_M1_EtoM1:
            return "MANNER_M1_EtoM1";
        case MANNER_M2toM1:
            return "MANNER_M2toM1";
        case MANNER_M2toMulti:
            return "MANNER_M2toMulti";
        case MANNER_S_M2toMulti:
            return "MANNER_S_M2toMulti";
        case MANNER_Multi_EtoMulti:
            return "MANNER_Multi_EtoMulti";
        case MANNER_MultitoC:
            return "MANNER_MultitoC";
        case MANNER_F_EtoF:
            return "MANNER_F_EtoF";
        case MANNER_CtoF:
            return "MANNER_CtoF";
        case MANNER_F_CtoF:
            return "MANNER_F_CtoF";
        case MANNER_CStoC:
            return "MANNER_CStoC";
        case MANNER_S_CStoC:
            return "MANNER_S_CStoC";
        case MANNER_C_StoCS:
            return "MANNER_C_StoCS";
        default:
            return "UNKNOWN";
    }
}

inline int GetIndex(int i, int nuci, int nucj){
    return i*square_base_number+nuci*base_number+nucj;
}

inline std::tuple<int, int, int> GetIndexTuple(int index_i_j) {
    int i = index_i_j / square_base_number;
    int nuc_index = index_i_j % square_base_number;
    int nuci = nuc_index / base_number;
    int nucj = nuc_index % base_number;

    return std::make_tuple(i, nuci, nucj);
}

struct State {
    int score;
    int index_1;
    int index_2;
    Manner MANNER;
    int last_pair_pos;

    State() : score(VALUE_MIN), index_1(-1), index_2(-1), MANNER(NONE), last_pair_pos(-1){}
    State(int score_, const int& index_1_, const int& index_2_, Manner MANNER_)
        : score(score_), index_1(index_1_), index_2(index_2_), MANNER(MANNER_){}
    State(int score_, const int& index_1_, const int& index_2_, Manner MANNER_, const int& last_pair_pos_)
        : score(score_), index_1(index_1_), index_2(index_2_), MANNER(MANNER_), last_pair_pos(last_pair_pos_){}
};

class AllTables {

public:

    std::vector<std::unordered_map<int, State>> bestN;
    std::vector<std::unordered_map<int, State>> bestS;
    std::vector<std::unordered_map<int, State>> bestF;
    std::vector<std::unordered_map<int, State>> bestC;
    std::vector<std::unordered_map<int, State>> bestCS;
    std::vector<std::unordered_map<int, State>> bestM1;
    std::vector<std::unordered_map<int, State>> bestM2;
    std::vector<std::unordered_map<int, State>> bestMulti;

    AllTables(std::string rna_seq, int size_) {
        bestN = initi_vec(rna_seq, size_);
        bestS = initi_vec(rna_seq, size_);
        bestF = initi_vec(rna_seq, size_);
        bestC = initi_vec(rna_seq, size_);
        bestCS = initi_vec(rna_seq, size_);
        bestM1 = initi_vec(rna_seq, size_);
        bestM2 = initi_vec(rna_seq, size_);
        bestMulti = initi_vec(rna_seq, size_);
    }

    std::vector<std::unordered_map<int, State>> initi_vec(std::string& rna_seq, int size_){
        std::vector<std::unordered_map<int, State>> seq_tables_(size_);
        return seq_tables_;
    }

};


int BASE(const std::string& x) {
    static const std::unordered_map<std::string, int> baseMap = {
        {"A", 1},
        {"C", 2},
        {"G", 3},
        {"U", 4},
        {"G_AG1", 5},
        {"G_AG2", 6},
        {"U_AG1", 7},
        {"U_AG2", 8},
        {"A1", 9},
        {"A2", 10},
        {"G1", 11},
        {"G2", 12},
        {"G_CU", 13},
        {"U1", 14},
        {"C1", 15}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return 0;
    }
}

static const int EnrBASEArr [] = { -1, 1, 2, 3, 4, 3, 3, 4, 4, 1, 1, 3, 3, 3, 4, 2 };

int EnrBASE(int x) {
    return EnrBASEArr[x];
    /* static const std::unordered_map<int, int> baseMap = {
        {1, 1}, //A
        {2, 2}, //C
        {3, 3}, //G
        {4, 4}, //U
        {5, 3}, //G_AG1
        {6, 3}, //G_AG2
        {7, 4}, //U_AG1
        {8, 4}, //U_AG2
        {9, 1}, //A1
        {10, 1}, //A2
        {11, 3}, //G1
        {12, 3}, //G2
        {13, 3}, //G_CU
        {14, 4}, //U1
        {15, 2} //C1
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return -1;
    } */
}


std::string CheckReBASE(int x) {
    static const std::unordered_map<int, std::string> baseMap = {
        {1, "A"},
        {2, "C"},
        {3, "G"},
        {4, "U"},
        {5, "G_AG1"},
        {6, "G_AG2"},
        {7, "U_AG1"},
        {8, "U_AG2"},
        {9, "A1"},
        {10, "A2"},
        {11, "G1"},
        {12, "G2"},
        {13, "G_CU"},
        {14, "U1"},
        {15, "C1"}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return "*";
    }
}

char reBASE(int x) {
    static const std::unordered_map<int, char> baseMap = {
        {1, 'A'},
        {2, 'C'},
        {3, 'G'},
        {4, 'U'},
        {5, 'G'},
        {6, 'G'},
        {7, 'U'},
        {8, 'U'},
        {9, 'A'},
        {10, 'A'},
        {11, 'G'},
        {12, 'G'},
        {13, 'G'},
        {14, 'U'},
        {15, 'C'},
        {-1, '*'}
    };

    auto it = baseMap.find(x);
    if (it != baseMap.end()) {
        return it->second;
    } else {
        return '*';
    }
}

void initialize(){
    // AU
    // UA
    std::vector<int> base_A_vec = {BASE("A"), BASE("A1"), BASE("A2")};
    std::vector<int> base_U_vec = {BASE("U"), BASE("U_AG1"), BASE("U_AG2"), BASE("U1")};

    for (int i = 0; i < base_A_vec.size(); i++){
        for (int j = 0; j < base_U_vec.size(); j++){
            _allowed_pairs[base_A_vec[i]][base_U_vec[j]] = true;
            _allowed_pairs[base_U_vec[j]][base_A_vec[i]] = true;
        }
    }

    //CG
    //GC
    std::vector<int> base_C_vec = {BASE("C"), BASE("C1")};
    std::vector<int> base_G_vec = {BASE("G"), BASE("G_AG1"), BASE("G_AG2"), BASE("G1"), BASE("G2"), BASE("G_CU")};
    
    for (int i = 0; i < base_C_vec.size(); i++){
        for (int j = 0; j < base_G_vec.size(); j++){
            _allowed_pairs[base_C_vec[i]][base_G_vec[j]] = true;
            _allowed_pairs[base_G_vec[j]][base_C_vec[i]] = true;
        }
    }

    //UG
    //GU
    for (int i = 0; i < base_U_vec.size(); i++){
        for (int j = 0; j < base_G_vec.size(); j++){
            _allowed_pairs[base_U_vec[i]][base_G_vec[j]] = true;
            _allowed_pairs[base_G_vec[j]][base_U_vec[i]] = true;
        }
    }
    
    // Arginine
    xcon_table[Arg][0][BASE("A")][BASE("G_AG1")] = true; //R
    xcon_table[Arg][0][BASE("C")][BASE("G_AG2")] = true; //R
    xcon_table[Arg][0][BASE("C")][BASE("G_CU")] = true; //R
    xcon_table[Arg][1][BASE("G_AG1")][BASE("A1")] = true; //R
    xcon_table[Arg][1][BASE("G_AG1")][BASE("G1")] = true; //R
    xcon_table[Arg][1][BASE("G_AG2")][BASE("A2")] = true; //R 
    xcon_table[Arg][1][BASE("G_AG2")][BASE("G2")] = true; //R
    xcon_table[Arg][1][BASE("G_CU")][BASE("C")] = true; //R
    xcon_table[Arg][1][BASE("G_CU")][BASE("U")] = true; //R

    // Leucine
    xcon_table[Leu][0][BASE("U")][BASE("U_AG1")] = true; //L
    xcon_table[Leu][0][BASE("C")][BASE("U_AG2")] = true; //L
    xcon_table[Leu][0][BASE("C")][BASE("U1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG1")][BASE("A1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG1")][BASE("G1")] = true; //L
    xcon_table[Leu][1][BASE("U_AG2")][BASE("A2")] = true; //L
    xcon_table[Leu][1][BASE("U_AG2")][BASE("G2")] = true; //L
    xcon_table[Leu][1][BASE("U1")][BASE("C")] = true; //L
    xcon_table[Leu][1][BASE("U1")][BASE("U")] = true; //L

    //Serine
    xcon_table[Ser][0][BASE("U")][BASE("C")] = true; //S
    xcon_table[Ser][0][BASE("A")][BASE("G")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("U")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("C")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("A")] = true; //S
    xcon_table[Ser][1][BASE("C")][BASE("G")] = true; //S
    xcon_table[Ser][1][BASE("G")][BASE("U1")] = true; //S
    xcon_table[Ser][1][BASE("G")][BASE("C1")] = true; //S

}

//...[i][j]...
bool IsLegal(std::vector<int>& con_seq, int base_i, int base_j, int pos_i){
    //j must be i+1
    int amino_ = con_seq[pos_i];
    int pos_j = pos_i+1;
    if((pos_i/3) == (pos_j/3)){
        if(amino_ != normal_ami){
            int amino_index = pos_i%3;
            if(xcon_table[amino_][amino_index][base_i][base_j]){
                return true;}else{return false;}
        }else{return true;}   
    }else{return true;}
    
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void ami_to_rna(std::vector<std::string>& rna_seq_list, std::string& ami_seq){
    std::string rna_seq;
    std::string Amino1;
    std::string codon;
    for(int i = 0; i < ami_seq.size(); i++){
        Amino1 = ami_seq[i];
        codon = Amino_label.at(Amino1);
        for (int j = 0; j < codon.size(); j++){
            rna_seq.push_back(codon[j]);
        }       
    }
    rna_seq_list.push_back(rna_seq);
}

void add_con_seq(std::vector<std::vector<int>>& con_seq_list, std::string& ami_seq){
    std::vector<int> con_seq;
    char Amino1; 
    for(int i=0; i<ami_seq.size(); i++){
        Amino1 = ami_seq[i];
        if(Amino1=='R'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Arg); //R==0
        }else if(Amino1=='L'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Leu); //L==1
        }else if(Amino1=='S'){
            for(int j=0; j<3; j++)
                con_seq.push_back(Ser); //S==2
        }else{
            for(int j=0; j<3; j++)
                con_seq.push_back(normal_ami); //-1
        }
    }
    con_seq_list.push_back(con_seq);
}

void check_valid_ami(int seq_size, std::string result_seq, std::string& ami_seq){
    std::string sub_str;
    std::string index;
    for(int i=0; i<seq_size; i++){
        sub_str =  result_seq.substr(3*i,3);
        if(re_Amino_label.count(sub_str)){
            index = re_Amino_label.at(sub_str);    
        }else{
            std::cout << "error aminoacid: ";
            std::cout << sub_str << std::endl;
            std::cout << "error index: ";
            std::cout << i << std::endl;
            break;
        }
                
        ami_seq.push_back(index[0]);
        
    }
}

inline void log_CAI(std::vector<float>& cai_vector, std::vector<float>& log_cai_vector){

    for(int i = 0; i < cai_vector.size(); i++){
        // log_cai_vector.push_back(log10f(cai_vector[i]));
        log_cai_vector.push_back(logf(cai_vector[i]));
    }
}

inline int getLastExtendedNuc(char amino_ ,std::string codon_){
    int last_nuc = -1;

    if(amino_ == 'R'){
        if(codon_ == "AGA")
            last_nuc = BASE("A1");
        if(codon_ == "AGG")
            last_nuc = BASE("G1");
        
        if(codon_ == "CGA")
            last_nuc = BASE("A2");
        if(codon_ == "CGG")
            last_nuc = BASE("G2");
        if(codon_ == "CGC")
            last_nuc = BASE("C");
        if(codon_ == "CGU")
            last_nuc = BASE("U");
        // {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
    }
    else if(amino_ == 'L'){
        if(codon_ == "CUU")
            last_nuc = BASE("U");
        if(codon_ == "CUC")
            last_nuc = BASE("C");
        if(codon_ == "CUA")
            last_nuc = BASE("A2");
        if(codon_ == "CUG")
            last_nuc = BASE("G2");
        
        if(codon_ == "UUA")
            last_nuc = BASE("A1");
        if(codon_ == "UUG")
            last_nuc = BASE("G1");
        // {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
    }
    else if(amino_ == 'S'){
        if(codon_ == "UCU")
            last_nuc = BASE("U");
        if(codon_ == "UCC")
            last_nuc = BASE("C");
        if(codon_ == "UCA")
            last_nuc = BASE("A");
        if(codon_ == "UCG")
            last_nuc = BASE("G");
        
        if(codon_ == "AGU")
            last_nuc = BASE("U1");
        if(codon_ == "AGC")
            last_nuc = BASE("C1");
        // {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
    }
    else{
        last_nuc = BASE(std::string(1, codon_[2]));
    }
    return last_nuc;
}

inline void initialize_CAI_table(std::vector<float> cai_vector_){
    std::vector<float> cai_vector;
    log_CAI(cai_vector_, cai_vector);


    for (int k = 0; k < CodonSet.size(); ++k) {
        char amino_ = CodonSet[k][0][0];

        if (CodonSetCAIMap.find(amino_) == CodonSetCAIMap.end()) {
            CodonSetCAIMap[amino_] = std::unordered_map<int, int>();
        }
    }


    for (int i = 0; i < CodonSet.size(); ++i){

        char amino_ = CodonSet[i][0][0];
        std::string codon_ = CodonSet[i][1];
        int last_nuc = BASE(std::string(1, codon_[2]));

        if(amino_ == 'R'){
            if(codon_ == "AGA")
                last_nuc = BASE("A1");
            if(codon_ == "AGG")
                last_nuc = BASE("G1");
            
            if(codon_ == "CGA")
                last_nuc = BASE("A2");
            if(codon_ == "CGG")
                last_nuc = BASE("G2");
            if(codon_ == "CGC")
                last_nuc = BASE("C");
            if(codon_ == "CGU")
                last_nuc = BASE("U");
            // {"CGU","R"}, {"CGC","R"}, {"CGA","R"}, {"CGG","R"}, {"AGA","R"}, {"AGG","R"},
        }
        else if(amino_ == 'L'){
            if(codon_ == "CUU")
                last_nuc = BASE("U");
            if(codon_ == "CUC")
                last_nuc = BASE("C");
            if(codon_ == "CUA")
                last_nuc = BASE("A2");
            if(codon_ == "CUG")
                last_nuc = BASE("G2");
            
            if(codon_ == "UUA")
                last_nuc = BASE("A1");
            if(codon_ == "UUG")
                last_nuc = BASE("G1");
            // {"CUU","L"}, {"CUC","L"}, {"CUA","L"}, {"CUG","L"}, {"UUA","L"}, {"UUG","L"},
        }
        else if(amino_ == 'S'){
            if(codon_ == "UCU")
                last_nuc = BASE("U");
            if(codon_ == "UCC")
                last_nuc = BASE("C");
            if(codon_ == "UCA")
                last_nuc = BASE("A");
            if(codon_ == "UCG")
                last_nuc = BASE("G");
            
            if(codon_ == "AGU")
                last_nuc = BASE("U1");
            if(codon_ == "AGC")
                last_nuc = BASE("C1");
            // {"UCU","S"}, {"UCC","S"}, {"UCA","S"}, {"UCG","S"}, {"AGU","S"}, {"AGC","S"},
        }
        float weight_ = LDCONST*lambda*cai_vector[i];
        CodonSetCAIMap[amino_][last_nuc] = int(weight_);
        
    }    

    
}

void check_ami_solution(std::string ami_seq, std::string rna_solution){
    std::string check_ami_seq;
    check_valid_ami(ami_seq.size(), rna_solution, check_ami_seq);


    // //check correct
    if(ami_seq != check_ami_seq){
        std::cout << "Current backtracked amino acids sequence:" << std::endl;
        std::cout << check_ami_seq << std::endl;
        std::cout << "ERROR: Codon sequences cannot transfer into amino acids." << std::endl;
    }
}