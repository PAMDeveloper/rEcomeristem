#ifndef UTILS_PARAMETERS_READER_HPP
#define UTILS_PARAMETERS_READER_HPP

#include <ModelParameters.hpp>

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


inline double round( double val, int decimal )
{
    double factor = std::pow(10, decimal);
    val *= factor;
    val = std::round(val);
    return val/factor;
}

using namespace ecomeristem;
using namespace std;
#include <iostream>
#include <fstream>
namespace utils {

class ParametersReader {
public:
    ParametersReader()
    { }

    virtual ~ParametersReader()
    { }

    void loadParametersFromFiles(const std::string &folder, ModelParameters &parameters) {
        std::ifstream varietyParams(folder + "/ECOMERISTEM_parameters.txt");
        std::string line;

        while (varietyParams >> line) {
            std::string s = line.substr(line.find("=") + 1, line.size());
            char* p;
            double converted = strtod(s.c_str(), &p);
            if (*p) {
                converted = JulianDayConverter::toJulianDayNumber(s, DATE_FORMAT_DMY);
            }
            parameters.set(line.substr(0, line.find("=")), converted);
        }

        std::ifstream * meteoFiles[5] = {
            new std::ifstream(folder + "/meteo_T.txt"),
            new std::ifstream(folder + "/meteo_PAR.txt"),
            new std::ifstream(folder + "/meteo_ETP.txt"),
            new std::ifstream(folder + "/meteo_irrig.txt"),
            new std::ifstream(folder + "/meteo_P.txt")
        };


        std::string values[5] = { "", "", "", "", "" };
        std::string date;
        bool first = true;

        while (*meteoFiles[0] >> date) {
            if (first) {
                double julianBegin = JulianDayConverter::toJulianDayNumber(date, DATE_FORMAT_DMY);
                parameters.set("BeginDate", julianBegin);
                parameters.beginDate = julianBegin;
                first = false;
            }

            *meteoFiles[0] >> values[0] >> values[0];

            for (int i = 1; i < 5; i++)
                *meteoFiles[i] >> values[i] >> values[i] >> values[i];


            //#ifdef OPTIM_NO_LEXCAST
            parameters.meteoValues.push_back(
                        Climate(
                            round(std::stod(values[0]),9),
                        round(std::stod(values[1]),9),
                    round(std::stod(values[2]),9),
                    round(std::stod(values[3]),9),
                    round(std::stod(values[4]),9)
                    )
                    );
            //#else
            //		    parameters.meteoValues.push_back(
            //		      model::models::Climate(
            //		        boost::lexical_cast<double>(values[0]),
            //		        boost::lexical_cast<double>(values[1]),
            //		        boost::lexical_cast<double>(values[2]),
            //		        boost::lexical_cast<double>(values[3]),
            //		        boost::lexical_cast<double>(values[4])
            //		      )
            //		    );
            //#endif
        }
        for (int i = 0; i < 5; i++)
          meteoFiles[i]->close();

        for (int i = 0; i < 5; ++i) {
            delete meteoFiles[i];
        }
        parameters.set("EndDate", JulianDayConverter::toJulianDayNumber(date, DATE_FORMAT_DMY));
        varietyParams.close();


    }

    map<string, vector<double>> loadVObsFromFile(const std::string &file_path) {
        std::ifstream vObsFile(file_path);
        std::string line;
        std::getline(vObsFile, line); //headers
        vector<string> headers;
        istringstream iss(line);
        copy(istream_iterator<string>(iss),
             istream_iterator<string>(),
             back_inserter(headers));

        map<string, vector<double> > obs;
        for (string h: headers) {
            string * s = new string(h);
            transform(s->begin(), s->end(), s->begin(), ::tolower);
            obs.insert ( std::pair<string,vector<double> >(*s, vector<double>()) );
            delete s;
        }

        while (std::getline(vObsFile, line))
        {
            vector<string> data = split(line, '\t');
            for (int i = 0; i < data.size(); ++i) {
                string s = data[i];
                string * h = new string(headers[i]);
                transform(h->begin(), h->end(), h->begin(), ::tolower);
                char* p;
                double converted = strtod(s.c_str(), &p);
                if (*p) {
                    obs[*h].push_back(nan(""));
                }
                else {
                    obs[*h].push_back(converted);
                }
                delete h;
            }
        }
        vObsFile.close();

        return obs;
    }


    //		void loadParametersFromProstgresql(const std::string &id,
    //            ModelParameters &parameters);

private:
    //        void load_meteo(PGconn* connection, ModelParameters &parameters);
    //		void load_data(PGconn* connection,
    //			const std::string &table,
    //			const std::string &id,
    //			const std::vector < std::string > &names,
    //            ModelParameters &parameters);

    //		void load_simulation(const std::string &id,
    //			PGconn* connection,
    //            ModelParameters &parameters);
};
} // namespace utils

#endif

