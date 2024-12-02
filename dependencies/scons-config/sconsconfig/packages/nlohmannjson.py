from .Package import Package

example = r"""
#include <iostream>
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

int main()
{
  // create a JSON object
  json j =
  {
    {"pi", 3.141},
    {"happy", true},
    {"name", "Niels"},
    {"nothing", nullptr},
    {
      "answer", {
        {"everything", 42}
      }
    },
    {"list", {1, 0, 2}},
    {
      "object", {
         {"currency", "USD"},
         {"value", 42.99}
      }
    }
  };

  // add new values
  j["new"]["key"]["value"] = {"another", "list"};

  // count elements
  auto s = j.size();
  j["size"] = s;

  // pretty print with indent of 4 spaces
  std::cout << std::setw(4) << j << '\n';
}
"""


class NlohmannJson(Package):
    def __init__(self, **kwargs):
        defaults = {
            "download_url": "https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.tar.gz"
        }
        defaults.update(kwargs)
        super(NlohmannJson, self).__init__(**defaults)
        self.sub_dirs = [
            ("", ""),
            ("include","")
        ]
        self.ext = '.cpp'

        self.check_text = example
        self.static = False
        self.number_output_lines = 11540 * 1.45

    def check(self, ctx):
        env = ctx.env
        ctx.Message("Checking for NlohmannJson ...  ")

        self.set_build_handler(
            [
                "mkdir -p ${PREFIX}",
                "cd ${SOURCE_DIR} && \
                 mkdir -pv build && cd build && \
                 cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=${PREFIX} .. && \
                 make -j && make install",
            ]
        )

        self.check_options(env)
        res = super(NlohmannJson, self).check(ctx)

        self.check_required(res[0], ctx)
        ctx.Result(res[0])
        return res[0]
