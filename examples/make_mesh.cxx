#include "NbtIO.hxx"
#include "NbsIO.hxx"
#include "MeshGenerator.hxx"

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static std::string trim(const std::string& s) {
    std::size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}

static void splitTokens(const std::string& line, std::vector<std::string>& out) {
    out.clear(); std::istringstream iss(line); std::string t; while (iss >> t) out.push_back(t);
}

static std::string dirnameOf(const std::string& path) {
    auto pos = path.find_last_of("/");
    if (pos == std::string::npos) return std::string(".");
    if (pos == 0) return std::string("/");
    return path.substr(0, pos);
}

static std::string joinPath(const std::string& a, const std::string& b) {
    if (b.empty()) return a;
    if (!b.empty() && (b[0] == '/'
#ifdef _WIN32
        || (b.size() > 1 && b[1] == ':')
#endif
        )) return b; // absolute
    if (a.empty() || a == ".") return b;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

static std::string replaceExt(const std::string& path, const std::string& newExt) {
    auto pos = path.find_last_of('.');
    if (pos == std::string::npos) return path + newExt;
    return path.substr(0, pos) + newExt;
}

// Read only the nbs_file entry from an .nbt file
static bool peekNbsPath(const std::string& nbtPath, std::string& outNbs) {
    std::ifstream ifs(nbtPath);
    if (!ifs) return false;
    std::string line; std::vector<std::string> toks;
    while (std::getline(ifs, line)) {
        auto hashPos = line.find('#'); if (hashPos != std::string::npos) line = line.substr(0, hashPos);
        std::string t = trim(line);
        if (t.empty() || t[0] == '*') continue;
        splitTokens(t, toks);
        if (!toks.empty() && toks[0] == "nbs_file" && toks.size() >= 2) {
            outNbs = toks[1];
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::fprintf(stderr, "Usage: %s <topology.nbt> <target_h> [out.msh]\n", argv[0]);
        return 2;
    }
    const std::string nbtPath = argv[1];
    const double h = std::atof(argv[2]);
    std::string mshPath = (argc >= 4) ? argv[3] : replaceExt(nbtPath, ".msh");

    // Get referenced .nbs file and resolve relative to nbt directory if needed
    std::string nbsRef;
    if (!peekNbsPath(nbtPath, nbsRef)) {
        std::fprintf(stderr, "Failed to read nbs_file from %s\n", nbtPath.c_str());
        return 1;
    }
    const std::string baseDir = dirnameOf(nbtPath);
    const std::string nbsPath = joinPath(baseDir, nbsRef);

    // Load curves
    std::vector<Nurbs> curves;
    std::string err;
    if (!NbsIO::readFile(nbsPath, curves, &err)) {
        std::fprintf(stderr, "Failed to read NBS %s: %s\n", nbsPath.c_str(), err.c_str());
        return 1;
    }

    // Load loops/regions using NbtIO (now that we have curves with ids)
    std::string nbsFromNbt;
    std::vector<Loop> loops; std::vector<Region> regions;
    if (!NbtIO::readFile(nbtPath, curves, nbsFromNbt, loops, regions, &err)) {
        std::fprintf(stderr, "Failed to read NBT %s: %s\n", nbtPath.c_str(), err.c_str());
        return 1;
    }

    // Generate mesh
    if (!MeshGenerator::generate(regions, curves, h, mshPath)) {
        std::fprintf(stderr, "Mesh generation failed\n");
        return 1;
    }
    std::printf("Wrote mesh: %s\n", mshPath.c_str());
    return 0;
}
