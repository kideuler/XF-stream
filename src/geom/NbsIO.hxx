#ifndef NBS_IO_HXX
#define NBS_IO_HXX

#include "Nurbs.hxx"
#include <string>
#include <vector>

namespace NbsIO {

// Read a .nbs file and append curves to 'out'. Returns true on success.
bool readFile(const std::string& path, std::vector<Nurbs>& out);
bool readFile(const std::string& path, std::vector<Nurbs>& out, std::string* errorMessage);

// Write curves to a .nbs file using explicit Nurbs representation (degree/knots/weights/control_points). Returns true on success.
bool writeFile(const std::string& path, const std::vector<Nurbs>& curves);
bool writeFile(const std::string& path, const std::vector<Nurbs>& curves, std::string* errorMessage);

} // namespace NbsIO

#endif // NBS_IO_HXX