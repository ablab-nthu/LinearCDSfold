#pragma once

#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <limits>
#include <cctype>
#include <chrono>
#include <ctime>
#include <math.h>

#include "utility.h"
#include "beamprune.h"


using namespace std;

template <typename T>
T BackTrack(string& rna_solution, string& structure_solution, string& rna_seq, AllTables<T>& alltables, vector<int>& con_seq){
    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::unordered_map<int, State<T>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    
    int seq_length = rna_seq.size();
    rna_solution.assign(seq_length, '.');
    structure_solution.assign(seq_length, '.');

    vector<int> first_nuclist = Base_table.at(rna_seq[0]);//i list
    vector<int> last_nuclist = Base_table.at(rna_seq[seq_length-1]);//j list

    State<T> state;
    T maxscore = VALUE_MIN<T>();
    int maxindex;
    for(int ci = 0; ci < first_nuclist.size(); ci++){
        for(int cj = 0; cj < last_nuclist.size(); cj++){
            int nuci = first_nuclist[ci];
            int nucj = last_nuclist[cj];
            int index = GetIndex(0, nuci, nucj);
            if(bestF[seq_length-1].count(index)){
                T newscore = bestF[seq_length-1][index].score;
                if (maxscore < newscore){
                    maxscore = newscore;
                    maxindex = index;
                    state = bestF[seq_length-1][index];}}
    }}

    stack<tuple<int, int, State<T>>> stk;
    stk.push(make_tuple(seq_length-1, maxindex, state));

    while ( !stk.empty() ) {

        tuple<int, int, State<T>> top = stk.top();
        int j = get<0>(top);
        int cindex = get<1>(top);
        state = get<2>(top);
        int index1 = state.index_1;
        int index2 = state.index_2;
        int i, nuci, nucj;
        std::tie(i, nuci, nucj) = GetIndexTuple(cindex);

        Manner MANNER = state.MANNER;
        State<T> next_state, next_state1, next_state2;
        stk.pop();
        int k, nuck, nucl;

        switch (MANNER){
            
            case NONE:
                break;
            case MANNER_NONEtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_NONEtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                break;
            case MANNER_N_EtoN:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_NtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestN[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_EtoS:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestS[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestC[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_S_CtoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestS[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                break;
            case MANNER_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state1 = bestCS[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state1));
                break;
            case MANNER_S_CStoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state2 = bestCS[j-1][index2];
                stk.push(make_tuple(j-1 ,index2, next_state2));
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestS[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                break;
            case MANNER_C_StoCS:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestC[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestS[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_CtoM1:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M1_CtoM2:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestM1[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_M1_EtoM1:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestM1[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_M2toM1:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_M2toMulti:
                next_state = bestM2[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_S_M2toMulti:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestS[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestM2[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            case MANNER_Multi_EtoMulti:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_MultitoC:
                rna_solution[i]=reBASE(nuci);
                rna_solution[j]=reBASE(nucj);
                structure_solution[i]='(';
                structure_solution[j]=')';
                next_state = bestMulti[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_F_EtoF:
                rna_solution[j]=reBASE(nucj);
                structure_solution[j]='.';
                next_state = bestF[j-1][index1];
                stk.push(make_tuple(j-1 ,index1, next_state));
                break;
            case MANNER_CtoF:
                next_state = bestC[j][index1];
                stk.push(make_tuple(j ,index1, next_state));
                break;
            case MANNER_F_CtoF:
                std::tie(k, nuck, nucl) = GetIndexTuple(index2);
                next_state1 = bestF[k-1][index1];
                stk.push(make_tuple(k-1 ,index1, next_state1));
                next_state2 = bestC[j][index2];
                stk.push(make_tuple(j ,index2, next_state2));
                break;
            default:
                cout << "Something error" << endl;
        }
    }

    return maxscore;
}


inline bool isLastNuc(int pos){
    if(pos%3 == 2)
        return true;
    return false;
}

inline int getCAI(char Amino_label, int nuc){
    return CodonSetCAIMap[Amino_label][nuc];
}

inline float getdernaCAI(char Amino_label, int nuc){
    return CodonSetCAIMap_derna[Amino_label][nuc];
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, T newscore, int preindex_1, Manner MANNER){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, -1, MANNER);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, T newscore, int preindex_1, int preindex_2, Manner MANNER){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, preindex_2, MANNER);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, T newscore, int preindex_1, Manner MANNER, int last_pair_pos){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, -1, MANNER, last_pair_pos);
        stateMap[index] = newState;
    }
}

template <typename T>
void update(std::unordered_map<int, State<T>>& stateMap, int index, T newscore, int preindex_1, int preindex_2, Manner MANNER, int last_pair_pos){

    auto it = stateMap.find(index);

    if (it == stateMap.end() || newscore > it->second.score) {
        State<T> newState(newscore, preindex_1, preindex_2, MANNER, last_pair_pos);
        stateMap[index] = newState;
    }
}

template <typename T>
void LCDSfoldCAI(AllTables<T>& alltables, std::string& rna_seq, std::vector<int>& con_seq, std::string& ami_seq){
    //LCDS
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::unordered_map<int, State<T>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0, -1, MANNER_NONEtoF);
        update(bestS[0], index, 0, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){
        
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::unordered_map<int, State<T>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        
        for (auto &nucj : nucj_list){
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && lambda)
                newscore += getCAI(ami_seq[j/3], nucj);
                
            update(bestS_j, index, newscore, -1, MANNER_NONEtoS);
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);

            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {
                int index_j_1 = itemN_j_1.first;
                State<T> stateN_j_1 = itemN_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateN_j_1.score;

                    if(isLastNuc(j) && lambda)
                        newscore += getCAI(ami_seq[j/3], nucj);
                    // N->N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    newscore = - v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score;
                                    if(isLastNuc(j) && lambda)
                                        newscore += getCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && lambda)
                                        newscore += getCAI(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }

        BeamPrune(con_seq, rna_seq, bestN_j, bestF, !lambda);

        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1];
            for (auto &itemS_j_1 : bestS_j_1) {
                int index_j_1 = itemS_j_1.first;
                State<T> stateS_j_1 = itemS_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) <= SINGLE_MAX_LEN){
                    newscore = stateS_j_1.score;
                    if(isLastNuc(j) && lambda)
                        newscore += getCAI(ami_seq[j/3], nucj);
                    // S->S
                    index = GetIndex(i, nuci, nucj);
                    update(bestS_j, index, newscore, index_j_1, MANNER_S_EtoS);
                }
            }
        }

        BeamPrune(con_seq, rna_seq, bestS_j, bestF, !lambda);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestCS_j_1 = bestCS[j_1];
            for (auto &itemCS_j_1 : bestCS_j_1) {

                int index_j_1 = itemCS_j_1.first;
                State<T> stateCS_j_1 = itemCS_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);

                int i = i1 - 1;
                if(i >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                    for(auto &nuci : nuci_list){
                        if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                            int pre_index_2 = stateCS_j_1.index_2;
                            int pre_index_1 = stateCS_j_1.index_1;
                            int s, nucs, nucl;
                            std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                            int s_1 = s - 1;
                            int a, nuca, nucs_1;
                            std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                            //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                            newscore = -v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score;
                            
                            if(isLastNuc(j) && lambda)
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && lambda)
                                newscore += getCAI(ami_seq[i/3], nuci);
                            
                            //CS->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CStoC);
                        }
                    }

                    std::unordered_map<int, State<T>>& bestS_i = bestS[i];
                    for (auto &itemS_i : bestS_i) {
                        int index_i = itemS_i.first;
                        State<T> stateS_i = itemS_i.second;
                        int k1, nuck1, nuci;
                        std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                        int k = k1 - 1;
                        if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                            vector<int> nuck_list = Base_table.at(rna_seq[k]); //i list
                            for(auto &nuck : nuck_list){
                                if(IsLegal(con_seq, nuck, nuck1, k) && _allowed_pairs[nuck][nucj]){
                                    int pre_index_2 = stateCS_j_1.index_2;
                                    int pre_index_1 = stateCS_j_1.index_1;
                                    int s, nucs, nucl;
                                    std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                                    int s_1 = s - 1;
                                    int a, nuca, nucs_1;
                                    std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);
                                    //v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);
                                    newscore = -v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                        EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)) + stateCS_j_1.score + stateS_i.score;
                                    
                                    if(isLastNuc(j) && lambda)
                                        newscore += getCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(k) && lambda)
                                        newscore += getCAI(ami_seq[k/3], nuck);
                                    
                                    //S+CS->C
                                    index = GetIndex(k, nuck, nucj);
                                    update(bestC_j, index, newscore, index_i ,index_j_1, MANNER_S_CStoC);
                                }
                            }
                        }

                    }
                }
            }
        }

        // BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemS_j : bestS_j) {
            int index_j = itemS_j.first;
            State<T> stateS_j = itemS_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestC_i_1 = bestC[i_1];
                
                for (auto &itemC_i_1 : bestC_i_1) {
                    int index_i_1 = itemC_i_1.first;
                    State<T> stateC_i_1 = itemC_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(k >= 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateS_j.score + stateC_i_1.score;
                        // C+S->CS
                        index = GetIndex(k, nuck, nucj);
                        update(bestCS_j, index, newscore, index_i_1, index_j, MANNER_C_StoCS);
                    }
                }
            }
        
        }

        BeamPrune(con_seq, rna_seq, bestCS_j, bestF, false);

        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                
                
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){
                            
                            newscore = -v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score;

                            if(isLastNuc(j) && lambda)
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && lambda)
                                newscore += getCAI(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        std::unordered_map<int, State<T>>& bestS_i = bestS[i];
                        for (auto &itemS_i : bestS_i) {
                            int index_i = itemS_i.first;
                            State<T> stateS_i = itemS_i.second;
                            int k1, nuck1, nuci;
                            std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                            int k = k1 - 1;
                            if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                for(auto &nuck : nuck_list){
                                    
                                    if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){
                                        
                                        newscore = -v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                        EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)) + stateC_j_1.score + stateS_i.score;
                                        
                                        if(isLastNuc(j) && lambda)
                                            newscore += getCAI(ami_seq[j/3], nucj);
                                        if(isLastNuc(k) && lambda)
                                            newscore += getCAI(ami_seq[k/3], nuck);
                                        
                                        //S+C->C
                                        index = GetIndex(k, nuck, nucj);
                                        update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                        }}}}
                    }
            }} //Cj_1
            
        } //nucj

        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) {
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && lambda)
                            newscore += getCAI(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            newscore = - v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)) + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && lambda)
                                newscore += getCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && lambda)
                                newscore += getCAI(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)) + stateC_j.score + stateM1_i_1.score;
                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && lambda)
                        newscore += getCAI(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1];
                for (auto &itemS_i_1 : bestS_i_1) {
                    int index_i_1 = itemS_i_1.first;
                    State<T> stateS_i_1 = itemS_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateM2_j.score + stateS_i_1.score;
                        //S+M2->Multi
                        index = GetIndex(k, nuck, nucj);
                        update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                    }
            }}
            
        }

        BeamPrune(con_seq, rna_seq, bestM1_j, bestF, false);
        BeamPrune(con_seq, rna_seq, bestMulti_j, bestF, false);

        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && lambda)
                        newscore += getCAI(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0){
                newscore = stateC_j.score - v_score_external_paired(nuci, nucj);
                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(k == 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateF_i_1.score + stateC_j.score - v_score_external_paired(nuci, nucj);
                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        BeamPrune(con_seq, rna_seq, bestF_j, bestF, false);
        }

    } //jend
}


template <typename T>
void LCDSfoldCAI_derna(AllTables<T>& alltables, string& rna_seq, vector<int>& con_seq, string& ami_seq){
    //LCDS
    int seq_length = rna_seq.size();
    T newscore;
    int index;
    int preindex;

    std::vector<std::unordered_map<int, State<T>>>& bestN = alltables.bestN;
    std::vector<std::unordered_map<int, State<T>>>& bestS = alltables.bestS;
    std::vector<std::unordered_map<int, State<T>>>& bestF = alltables.bestF;
    std::vector<std::unordered_map<int, State<T>>>& bestC = alltables.bestC;
    std::vector<std::unordered_map<int, State<T>>>& bestCS = alltables.bestCS;
    std::vector<std::unordered_map<int, State<T>>>& bestM1 = alltables.bestM1;
    std::vector<std::unordered_map<int, State<T>>>& bestM2 = alltables.bestM2;
    std::vector<std::unordered_map<int, State<T>>>& bestMulti = alltables.bestMulti;

    vector<int> nuc0_list = Base_table.at(rna_seq[0]); //j list
    for(int c0 = 0; c0 < nuc0_list.size(); ++c0){
        int nuc0 = nuc0_list[c0];
        index = GetIndex(0, nuc0, nuc0);
        update(bestN[0], index, 0.0f, -1, MANNER_NONEtoN);
        update(bestF[0], index, 0.0f, -1, MANNER_NONEtoF);
        update(bestS[0], index, 0.0f, -1, MANNER_NONEtoS);
    }
    
    for(int j = 1; j < seq_length; ++j){
        std::unordered_map<int, State<T>>& bestN_j = bestN[j];
        std::unordered_map<int, State<T>>& bestS_j = bestS[j];
        std::unordered_map<int, State<T>>& bestF_j = bestF[j];
        std::unordered_map<int, State<T>>& bestC_j = bestC[j];
        std::unordered_map<int, State<T>>& bestCS_j = bestCS[j];
        std::unordered_map<int, State<T>>& bestM1_j = bestM1[j];
        std::unordered_map<int, State<T>>& bestM2_j = bestM2[j];
        std::unordered_map<int, State<T>>& bestMulti_j = bestMulti[j];
        
        vector<int> nucj_list = Base_table.at(rna_seq[j]); //j list
        int j_1 = j - 1;
        
        for (auto &nucj : nucj_list){
            
            // NONE->S
            // NONE->N
            index = GetIndex(j, nucj, nucj);
            newscore = 0;

            if(isLastNuc(j) && lambda)
                newscore += getdernaCAI(ami_seq[j/3], nucj);
                
            update(bestS_j, index, newscore, -1, MANNER_NONEtoS);
            update(bestN_j, index, newscore, -1, MANNER_NONEtoN);

            std::unordered_map<int, State<T>>& bestN_j_1 = bestN[j_1];
            for (auto &itemN_j_1 : bestN_j_1) {
                int index_j_1 = itemN_j_1.first;
                State<T> stateN_j_1 = itemN_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateN_j_1.score;

                    if(isLastNuc(j) && lambda)
                        newscore += getdernaCAI(ami_seq[j/3], nucj);
                    // N->N
                    index = GetIndex(i, nuci, nucj);
                    update(bestN_j, index, newscore, index_j_1, MANNER_N_EtoN);

                    int i_1 = i - 1;
                    if((i_1 >= 0) && (j - i_1 > HAIRPIN_GAP)){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                            if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                                if(_allowed_pairs[nuci_1][nucj]){
                                    // newscore = - v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)) + stateN_j_1.score; // for MFE
                                    int mfe_ = - v_score_hairpin(i_1, j, EnrBASE(nuci_1), EnrBASE(nuci), EnrBASE(nucj_1), EnrBASE(nucj)); // for MFE
                                    newscore = (lambda) * mfe_ + stateN_j_1.score;

                                    if(isLastNuc(j) && lambda)
                                        newscore += getdernaCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(i_1) && lambda)
                                        newscore += getdernaCAI(ami_seq[i_1/3], nuci_1);
                                    //N->C
                                    index = GetIndex(i_1, nuci_1, nucj);
                                    update(bestC_j, index, newscore, index_j_1, MANNER_NtoC);
                            }}
                }}}
            }
        }

        BeamPrune(con_seq, rna_seq, bestN_j, bestF, !lambda);

        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestS_j_1 = bestS[j_1];
            for (auto &itemS_j_1 : bestS_j_1) {
                int index_j_1 = itemS_j_1.first;
                State<T> stateS_j_1 = itemS_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (j - i) <= SINGLE_MAX_LEN){
                    newscore = stateS_j_1.score;
                    if(isLastNuc(j) && lambda)
                        newscore += getdernaCAI(ami_seq[j/3], nucj);
                    // S->S
                    index = GetIndex(i, nuci, nucj);
                    update(bestS_j, index, newscore, index_j_1, MANNER_S_EtoS);
                }
            }
        }

        BeamPrune(con_seq, rna_seq, bestS_j, bestF, !lambda);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestCS_j_1 = bestCS[j_1];
            for (auto &itemCS_j_1 : bestCS_j_1) {

                int index_j_1 = itemCS_j_1.first;
                State<T> stateCS_j_1 = itemCS_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);

                int i = i1 - 1;
                if(i >= 0 && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i list
                    for(auto &nuci : nuci_list){
                        if(IsLegal(con_seq, nuci, nuci1, i) && _allowed_pairs[nuci][nucj]){
                            int pre_index_2 = stateCS_j_1.index_2;
                            int pre_index_1 = stateCS_j_1.index_1;
                            int s, nucs, nucl;
                            std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                            int s_1 = s - 1;
                            int a, nuca, nucs_1;
                            std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);

                            int mfe_ = -v_score_single(i, j, i1, s_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                                EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)); // for MFE
                            newscore = (lambda) * mfe_ + stateCS_j_1.score;

                            if(isLastNuc(j) && lambda)
                                newscore += getdernaCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && lambda)
                                newscore += getdernaCAI(ami_seq[i/3], nuci);
                            
                            //CS->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CStoC);
                        }
                    }

                    std::unordered_map<int, State<T>>& bestS_i = bestS[i];
                    for (auto &itemS_i : bestS_i) {
                        int index_i = itemS_i.first;
                        State<T> stateS_i = itemS_i.second;
                        int k1, nuck1, nuci;
                        std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                        int k = k1 - 1;
                        if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                            vector<int> nuck_list = Base_table.at(rna_seq[k]); //i list
                            for(auto &nuck : nuck_list){
                                if(IsLegal(con_seq, nuck, nuck1, k) && _allowed_pairs[nuck][nucj]){
                                    int pre_index_2 = stateCS_j_1.index_2;
                                    int pre_index_1 = stateCS_j_1.index_1;
                                    int s, nucs, nucl;
                                    std::tie(s, nucs, nucl) = GetIndexTuple(pre_index_2);
                                    int s_1 = s - 1;
                                    int a, nuca, nucs_1;
                                    std::tie(a, nuca, nucs_1) = GetIndexTuple(pre_index_1);

                                    int mfe_ = -v_score_single(k, j, i1, s_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                    EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucs_1), EnrBASE(nucs)); // for MFE
                                    newscore = (lambda) * mfe_ + stateCS_j_1.score + stateS_i.score;

                                    if(isLastNuc(j) && lambda)
                                        newscore += getdernaCAI(ami_seq[j/3], nucj);
                                    if(isLastNuc(k) && lambda)
                                        newscore += getdernaCAI(ami_seq[k/3], nuck);
                                    
                                    //S+CS->C
                                    index = GetIndex(k, nuck, nucj);
                                    update(bestC_j, index, newscore, index_i ,index_j_1, MANNER_S_CStoC);
                                }
                            }
                        }

                    }
                }
            }
        }

        for (auto &itemS_j : bestS_j) {
            int index_j = itemS_j.first;
            State<T> stateS_j = itemS_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestC_i_1 = bestC[i_1];
                
                for (auto &itemC_i_1 : bestC_i_1) {
                    int index_i_1 = itemC_i_1.first;
                    State<T> stateC_i_1 = itemC_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(k >= 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateS_j.score + stateC_i_1.score;
                        // C+S->CS
                        index = GetIndex(k, nuck, nucj);
                        update(bestCS_j, index, newscore, index_i_1, index_j, MANNER_C_StoCS);
                    }
                }
            }
        
        }

        BeamPrune(con_seq, rna_seq, bestCS_j, bestF, false);

        for (auto &nucj : nucj_list){

            std::unordered_map<int, State<T>>& bestC_j_1 = bestC[j_1];
            for (auto &itemC_j_1 : bestC_j_1) {
                int index_j_1 = itemC_j_1.first;
                State<T> stateC_j_1 = itemC_j_1.second;
                int i1, nucj_1, nuci1;
                std::tie(i1, nuci1, nucj_1) = GetIndexTuple(index_j_1);
                int i = i1 - 1;

                
                
                if(IsLegal(con_seq, nucj_1, nucj, j_1) && (i >= 0) && (j - i > (HAIRPIN_GAP + 2))){
                    
                    vector<int> nuci_list = Base_table.at(rna_seq[i]); //i-1 list
                    for(auto &nuci : nuci_list){
                        
                        if(_allowed_pairs[nuci][nucj] && IsLegal(con_seq, nuci, nuci1, i) && IsLegal(con_seq, nucj_1, nucj, j_1)){

                            int mfe_ = -v_score_single(i, j, i1, j_1, EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj),
                            EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)); // for MFE

                            newscore = (lambda) * mfe_ + stateC_j_1.score;

                            if(isLastNuc(j) && lambda)
                                newscore += getdernaCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i) && lambda)
                                newscore += getdernaCAI(ami_seq[i/3], nuci);
                            
                            //C->C
                            index = GetIndex(i, nuci, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_CtoC);
                            
                            
                    }}

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        std::unordered_map<int, State<T>>& bestS_i = bestS[i];
                        for (auto &itemS_i : bestS_i) {
                            int index_i = itemS_i.first;
                            State<T> stateS_i = itemS_i.second;
                            int k1, nuck1, nuci;
                            std::tie(k1, nuck1, nuci) = GetIndexTuple(index_i);
                            int k = k1 - 1;
                            if(k >= 0 && IsLegal(con_seq, nuci, nuci1, i)){
                                vector<int> nuck_list = Base_table.at(rna_seq[k]); //k list
                                for(auto &nuck : nuck_list){
                                    
                                    if(_allowed_pairs[nuck][nucj] && IsLegal(con_seq, nuck, nuck1, k)){

                                        int mfe_ = -v_score_single(k, j, i1, j_1, EnrBASE(nuck), EnrBASE(nuck1), EnrBASE(nucj_1), EnrBASE(nucj),
                                        EnrBASE(nuci), EnrBASE(nuci1), EnrBASE(nucj_1), EnrBASE(nucj)); // for MFE
                                        newscore = (lambda) * mfe_ + stateC_j_1.score + stateS_i.score;

                                        if(isLastNuc(j) && lambda)
                                            newscore += getdernaCAI(ami_seq[j/3], nucj);
                                        if(isLastNuc(k) && lambda)
                                            newscore += getdernaCAI(ami_seq[k/3], nuck);
                                        
                                        //S+C->C
                                        index = GetIndex(k, nuck, nucj);
                                        update(bestC_j, index, newscore, index_i, index_j_1, MANNER_S_CtoC);

                        }}}}
                    }
            }} //Cj_1
            
        } //nucj

        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestMulti_j_1 = bestMulti[j_1];
            for (auto &itemMulti_j_1 : bestMulti_j_1) {
                int index_j_1 = itemMulti_j_1.first;
                State<T> stateMulti_j_1 = itemMulti_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    
                    int last_pair_pos_Multi = stateMulti_j_1.last_pair_pos;
                    if(j - last_pair_pos_Multi <= SINGLE_MAX_LEN){
                        
                        newscore = stateMulti_j_1.score;

                        if(isLastNuc(j) && lambda)
                            newscore += getdernaCAI(ami_seq[j/3], nucj);

                        //Multi->Multi
                        index = GetIndex(i, nuci, nucj);
                        update(bestMulti_j, index, newscore, index_j_1, MANNER_Multi_EtoMulti, last_pair_pos_Multi);
                    }

                    int i_1 = i - 1;
                    if(i_1 >= 0){
                        vector<int> nuci_1_list = Base_table.at(rna_seq[i_1]); //i-1 list
                        for(auto &nuci_1 : nuci_1_list){
                        if(IsLegal(con_seq, nuci_1, nuci, i_1) && _allowed_pairs[nuci_1][nucj]){
                            int mfe_ = - v_score_multi(EnrBASE(nuci_1), EnrBASE(nucj)); // for MFE

                            newscore = (lambda) * mfe_ + stateMulti_j_1.score;
                            
                            if(isLastNuc(j) && lambda)
                                newscore += getdernaCAI(ami_seq[j/3], nucj);
                            if(isLastNuc(i_1) && lambda)
                                newscore += getdernaCAI(ami_seq[i_1/3], nuci_1);

                            //Multi->C
                            index = GetIndex(i_1, nuci_1, nucj);
                            update(bestC_j, index, newscore, index_j_1, MANNER_MultitoC);
                        }}
                    }
                }
            }
        } //nucj
        
        BeamPrune(con_seq, rna_seq, bestC_j, bestF, false);

        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            int mfe_ = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)); // for MFE
            newscore = (lambda) * mfe_ + stateC_j.score;

            //C->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_CtoM1, j);

            int i_1 = i - 1;
            if(i_1 >= HAIRPIN_GAP + 1){
                std::unordered_map<int, State<T>>& bestM1_i_1 = bestM1[i_1];
                for (auto &itemM1_i_1 : bestM1_i_1) {
                    int index_i_1 = itemM1_i_1.first;
                    State<T> stateM1_i_1 = itemM1_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        int mfe_ = - v_score_M1(EnrBASE(nuci), EnrBASE(nucj)); // for MFE
                        newscore = (lambda) * mfe_ + stateC_j.score + stateM1_i_1.score;

                        // M1+C->M2
                        index = GetIndex(k, nuck, nucj);
                        update(bestM2_j, index, newscore, index_i_1 ,index_j, MANNER_M1_CtoM2, j);
                    }
                }
            }
        }
        
        BeamPrune(con_seq, rna_seq, bestM2_j, bestF, false);

        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestM1_j_1 = bestM1[j_1];
            for (auto &itemM1_j_1 : bestM1_j_1) {
                int index_j_1 = itemM1_j_1.first;
                State<T> stateM1_j_1 = itemM1_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);

                int last_pair_pos_M1 = stateM1_j_1.last_pair_pos;
                assert(last_pair_pos_M1 != -1);
                if(j - last_pair_pos_M1 <= SINGLE_MAX_LEN && IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateM1_j_1.score;

                    if(isLastNuc(j) && lambda)
                        newscore += getdernaCAI(ami_seq[j/3], nucj);

                    //M1->M1
                    index = GetIndex(i, nuci, nucj);
                    update(bestM1_j, index, newscore, index_j_1, MANNER_M1_EtoM1, last_pair_pos_M1);
                }
            }
        } //nucj

        for (auto &itemM2_j : bestM2_j) {
            int index_j = itemM2_j.first;
            State<T> stateM2_j = itemM2_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);
            newscore = stateM2_j.score;
            int last_pair_pos_M2 = stateM2_j.last_pair_pos;

            //M2->M1
            update(bestM1_j, index_j, newscore, index_j, MANNER_M2toM1, last_pair_pos_M2);
            //M2->Multi
            update(bestMulti_j, index_j, newscore, index_j, MANNER_M2toMulti, last_pair_pos_M2);
            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestS_i_1 = bestS[i_1];
                for (auto &itemS_i_1 : bestS_i_1) {
                    int index_i_1 = itemS_i_1.first;
                    State<T> stateS_i_1 = itemS_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(IsLegal(con_seq, nuci_1, nuci, i_1)){
                        newscore = stateM2_j.score + stateS_i_1.score;
                        //S+M2->Multi
                        index = GetIndex(k, nuck, nucj);
                        update(bestMulti_j, index, newscore, index_i_1, index_j, MANNER_S_M2toMulti, last_pair_pos_M2);
                    }
            }}
            
        }

        BeamPrune(con_seq, rna_seq, bestM1_j, bestF, false);
        BeamPrune(con_seq, rna_seq, bestMulti_j, bestF, false);

        {//conclusion!
        
        for (auto &nucj : nucj_list){
            std::unordered_map<int, State<T>>& bestF_j_1 = bestF[j_1];
            for (auto &itemF_j_1 : bestF_j_1) {
                int index_j_1 = itemF_j_1.first;
                State<T> stateF_j_1 = itemF_j_1.second;
                int i, nucj_1, nuci;
                std::tie(i, nuci, nucj_1) = GetIndexTuple(index_j_1);
                if(IsLegal(con_seq, nucj_1, nucj, j_1)){
                    newscore = stateF_j_1.score;
                    if(isLastNuc(j) && lambda)
                        newscore += getdernaCAI(ami_seq[j/3], nucj);
                    //F->F
                    index = GetIndex(i, nuci, nucj);
                    update(bestF_j, index, newscore, index_j_1, MANNER_F_EtoF);
                }
            }
        }
        
        
        for (auto &itemC_j : bestC_j) {
            int index_j = itemC_j.first;
            State<T> stateC_j = itemC_j.second;
            int i, nucj, nuci;
            std::tie(i, nuci, nucj) = GetIndexTuple(index_j);

            if(i == 0){
                int mfe_ = - v_score_external_paired(nuci, nucj); // for MFE
                newscore = (lambda) * mfe_ + stateC_j.score;

                //C->F
                update(bestF_j, index_j, newscore, index_j, MANNER_CtoF);
            }

            int i_1 = i - 1;
            if(i_1 >= 0){
                std::unordered_map<int, State<T>>& bestF_i_1 = bestF[i_1];
                for (auto &itemF_i_1 : bestF_i_1) {
                    int index_i_1 = itemF_i_1.first;
                    State<T> stateF_i_1 = itemF_i_1.second;
                    int k, nuck, nuci_1;
                    std::tie(k, nuck, nuci_1) = GetIndexTuple(index_i_1);
                    if(k == 0 && IsLegal(con_seq, nuci_1, nuci, i_1)){
                        int mfe_ = - v_score_external_paired(nuci, nucj); // for MFE
                        newscore = (lambda) * mfe_ + stateF_i_1.score + stateC_j.score;

                        //F+C->F
                        index = GetIndex(k, nuck, nucj);
                        update(bestF_j, index, newscore, index_i_1, index_j, MANNER_F_CtoF);
            }}}
        }
        
        BeamPrune(con_seq, rna_seq, bestF_j, bestF, false);
        }

    } //jend
}