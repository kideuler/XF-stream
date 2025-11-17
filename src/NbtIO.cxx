#include "NbtIO.hxx"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>

namespace {
static inline std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}
static void splitTokens(const std::string& line, std::vector<std::string>& out) {
    out.clear(); std::istringstream iss(line); std::string t; while (iss >> t) out.push_back(t);
}
}

namespace NbtIO {

bool writeFile(const std::string& path,
               const std::string& nbsFileName,
               const std::vector<Nurbs>& curves,
               const std::vector<Loop>& loops,
               const std::vector<Region>& regions,
               std::string* errorMessage) {
    try {
        std::ofstream ofs(path);
        if (!ofs) throw std::runtime_error("Could not open nbt file for writing");
        ofs << "* XF-stream NBT topology v1\n";
        ofs << "nbs_file " << nbsFileName << "\n";
        ofs << "loops " << loops.size() << "\n";
        for (std::size_t i = 0; i < loops.size(); ++i) {
            ofs << "loop " << i << " curves";
            for (const auto& oc : loops[i]) {
                if (oc.index >= curves.size()) throw std::runtime_error("Loop references curve index out of range during write");
                int cid = curves[oc.index].id();
                if (cid < 0) throw std::runtime_error("Curve missing id in topology write");
                int signedCid = oc.reversed ? -cid : cid;
                ofs << ' ' << signedCid;
            }
            ofs << '\n';
        }
        ofs << "regions " << regions.size() << "\n";
        for (std::size_t i = 0; i < regions.size(); ++i) {
            ofs << "region " << i << " outer";
            auto loopsEqual = [](const Loop& A, const Loop& B) {
                if (A.size() != B.size()) return false;
                for (std::size_t k = 0; k < A.size(); ++k) {
                    if (A[k].index != B[k].index || A[k].reversed != B[k].reversed) return false;
                }
                return true;
            };
            std::size_t outerIdx = loops.size();
            for (std::size_t L = 0; L < loops.size(); ++L) {
                if (loopsEqual(loops[L], regions[i].outer)) { outerIdx = L; break; }
            }
            if (outerIdx == loops.size()) throw std::runtime_error("Could not match outer loop content in region write");
            ofs << ' ' << outerIdx << " inners";
            for (const auto& inner : regions[i].inners) {
                std::size_t innerIdx = loops.size();
                for (std::size_t L = 0; L < loops.size(); ++L) {
                    if (loopsEqual(loops[L], inner)) { innerIdx = L; break; }
                }
                if (innerIdx == loops.size()) throw std::runtime_error("Could not match inner loop content in region write");
                ofs << ' ' << innerIdx;
            }
            ofs << '\n';
        }
        ofs << "end\n";
        return true;
    } catch (const std::exception& e) {
        if (errorMessage) *errorMessage = e.what();
        return false;
    }
}

bool readFile(const std::string& path,
              const std::vector<Nurbs>& curves,
              std::string& outNbsFileName,
              std::vector<Loop>& outLoops,
              std::vector<Region>& outRegions,
              std::string* errorMessage) {
    try {
        std::ifstream ifs(path);
        if (!ifs) throw std::runtime_error("Could not open nbt file for reading");
        std::string line; std::vector<std::string> toks;
        enum class State { Idle, InBody };
        State st = State::Idle;
        outLoops.clear(); outRegions.clear(); outNbsFileName.clear();
        while (std::getline(ifs, line)) {
            auto hashPos = line.find('#'); if (hashPos != std::string::npos) line = line.substr(0, hashPos);
            std::string t = trim(line); if (t.empty() || t[0] == '*') continue;
            splitTokens(t, toks); if (toks.empty()) continue;
            if (toks[0] == "nbs_file" && toks.size() >= 2) {
                outNbsFileName = toks[1];
            } else if (toks[0] == "loops" && toks.size() >= 2) {
                // reserve loops count but loops defined individually
            } else if (toks[0] == "loop") {
                // loop <i> curves <list>
                Loop L;
                for (std::size_t i = 0; i < toks.size(); ++i) {
                    if (toks[i] == "curves") {
                        for (std::size_t j = i+1; j < toks.size(); ++j) {
                            int signedCid = std::stoi(toks[j]);
                            bool reversed = signedCid < 0;
                            int cid = reversed ? -signedCid : signedCid;
                            // find curve with this id
                            std::size_t idx = curves.size();
                            for (std::size_t k = 0; k < curves.size(); ++k) {
                                if (curves[k].id() == cid) { idx = k; break; }
                            }
                            if (idx == curves.size()) throw std::runtime_error("Loop references unknown curve id");
                            L.push_back(OrientedCurve{idx, reversed});
                        }
                        break;
                    }
                }
                outLoops.push_back(std::move(L));
            } else if (toks[0] == "regions" && toks.size() >= 2) {
                // ignore count
            } else if (toks[0] == "region") {
                // region <i> outer <oidx> inners <i0> <i1> ...
                std::size_t outerLoopIdx = static_cast<std::size_t>(-1);
                std::vector<std::size_t> innerIdxs;
                for (std::size_t i = 0; i < toks.size(); ++i) {
                    if (toks[i] == "outer" && i+1 < toks.size()) {
                        outerLoopIdx = static_cast<std::size_t>(std::stol(toks[i+1]));
                    } else if (toks[i] == "inners") {
                        for (std::size_t j = i+1; j < toks.size(); ++j) {
                            innerIdxs.push_back(static_cast<std::size_t>(std::stol(toks[j])));
                        }
                        break;
                    }
                }
                if (outerLoopIdx == static_cast<std::size_t>(-1) || outerLoopIdx >= outLoops.size()) {
                    throw std::runtime_error("Region outer loop index invalid");
                }
                Region R; R.outer = outLoops[outerLoopIdx];
                for (auto li : innerIdxs) {
                    if (li >= outLoops.size()) throw std::runtime_error("Region inner loop index invalid");
                    R.inners.push_back(outLoops[li]);
                }
                outRegions.push_back(std::move(R));
            } else if (toks[0] == "end") {
                break;
            }
        }
        return true;
    } catch (const std::exception& e) {
        if (errorMessage) *errorMessage = e.what();
        return false;
    }
}

} // namespace NbtIO
