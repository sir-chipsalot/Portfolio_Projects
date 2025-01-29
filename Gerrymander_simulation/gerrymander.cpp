#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>
#include <tuple>

using namespace std;

// Voter class: Represents a voter with an ID and party affiliation
class Voter {
private:
    int id;  // Unique ID for the voter
    int affiliation; // 1 for Party A, -1 for Party B, 0 for undecided

public:
    Voter(int voter_id, int party_affiliation) : id(voter_id), affiliation(party_affiliation) {}
    int get_affiliation() const { return affiliation; }
};

// District class: Represents a district containing voters
class District {
private:
    vector<Voter> voters;

public:
    District(const vector<Voter>& voter_list) : voters(voter_list) {}

    vector<Voter>& get_voters() { return voters; }

    int majority() const {
        int countA = 0, countB = 0;

        for (const Voter& voter : voters) {
            if (voter.get_affiliation() == 1) countA++;
            else if (voter.get_affiliation() == -1) countB++;
        }

        if (countA > countB) return 1;  // Party A wins
        if (countB > countA) return -1; // Party B wins
        return 0;                       // Tie
    }

    string print() const {
        int majority_affiliation = majority();
        if (majority_affiliation == 1) return "Party A";
        if (majority_affiliation == -1) return "Party B";
        return "Tie";
    }
};

// Districting class: Manages the districting process
class Districting {
private:
    vector<District> districts;

public:
    Districting(vector<Voter> population, int district_size) {
        int num_districts = population.size() / district_size;

        for (int i = 0; i < num_districts; ++i) {
            vector<Voter> district_voters;
            for (int j = i * district_size; j < (i + 1) * district_size; ++j) {
                district_voters.push_back(population[j]);
            }
            districts.push_back(District(district_voters));
        }
    }

    void print() const {
        for (size_t i = 0; i < districts.size(); ++i) {
            cout << "District " << i + 1 << ": " << districts[i].print() << endl;
        }
    }

    tuple<int, int, int> lean() const {
        int countA = 0, countB = 0, countTies = 0;

        for (const District& district : districts) {
            int majority = district.majority();
            if (majority == 1) countA++;
            else if (majority == -1) countB++;
            else countTies++;
        }

        return make_tuple(countA, countB, countTies);
    }

    void gerrymander() {
        // Count the number of districts won by Party A and Party B
        int countA = 0, countB = 0;

        for (int i = 0; i < districts.size(); ++i) {
            int majority_affiliation = districts[i].majority();
            if (majority_affiliation == 1) {
                countA++;
            } else if (majority_affiliation == -1) {
                countB++;
            }
        }

        // Determine which party has fewer districts (the losing party)
        int losing_party = 0;
        if (countA < countB) {
            losing_party = 1; // Party A is losing
        } else if (countB < countA) {
            losing_party = -1; // Party B is losing
        }

        // Move voters from the winning party to the next district in order to manipulate the outcome
        for (int i = 0; i < districts.size() - 1; ++i) {
            vector<Voter>& voters_i = districts[i].get_voters();

            // If the current district is won by the winning party, try to move voters from it
            if (districts[i].majority() != losing_party && i + 1 < districts.size()) {
                vector<Voter>& voters_next = districts[i + 1].get_voters();

                int move_count = voters_i.size() / 4;  // Decide to move 25% of the voters in the district
                int moved = 0;

                // Move voters from the current district to the next district
                for (int j = 0; j < voters_i.size() && moved < move_count; ++j) {
                    if (voters_i[j].get_affiliation() == losing_party) {
                        voters_next.push_back(voters_i[j]); // Move the voter to the next district
                        voters_i.erase(voters_i.begin() + j);  // Remove the voter from the current district
                        moved++;
                        j--; // Adjust the index after removing a voter
                    }
                }
            }
        }
    }
};

int main() {
    srand(time(0));

    int population_size = 10000000;
    int district_size = 200000;
    vector<Voter> population;

    // Generate voters
    for (int i = 0; i < population_size; ++i) {
        int affiliation;
        int roll = rand() % 100;

        if (roll < 45) affiliation = 1;       // 45% chance for Party A
        else if (roll < 90) affiliation = -1; // 45% chance for Party B
        else affiliation = 0;                // 10% chance for undecided

        population.push_back(Voter(i, affiliation));
    }

    // Create the original (regular) districting
    cout << "Regular Districting Results:\n";
    Districting regular_districting(population, district_size);
    regular_districting.print();

    // Print summary of regular districting
    auto [regA, regB, regTies] = regular_districting.lean();
    cout << "\nSummary of Regular Districting:\n";
    cout << "Party A: " << regA << " districts, Party B: " << regB << " districts, Ties: " << regTies << endl;

    // Gerrymander the existing districts
    cout << "\nGerrymandered Districting Results:\n";
    regular_districting.gerrymander();
    regular_districting.print();

    // Print summary of gerrymandered districting
    auto [gerrA, gerrB, gerrTies] = regular_districting.lean();
    cout << "\nSummary After Gerrymandering:\n";
    cout << "Party A: " << gerrA << " districts, Party B: " << gerrB << " districts, Ties: " << gerrTies << endl;

    return 0;
}

