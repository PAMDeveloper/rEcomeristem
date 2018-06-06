#include <Rcpp.h>
#include <Rinternals.h>

#include "recomeristem_types.hpp"

using namespace Rcpp;


DataFrame mapOfVectorToDF(std::map<std::string, std::vector<double>> map) {
  Rcpp::CharacterVector names;
  List values(map.size());
  int idx = 0;
  for(auto const& it: map){
    names.push_back(it.first);
    NumericVector vec( it.second.begin(), it.second.end() );
    values[idx] = vec;
    idx++;
  }
  DataFrame df(values);
  df.attr("names") = names;
  return df;
}

map <string, vector <double > > mapFromDF(DataFrame list) {
  map <string, vector <double > > map;
  CharacterVector names = list.attr("names");
  for (int i = 0; i < names.size(); ++i) {
    std::vector<double> vec = as<std::vector<double> >(list[i]);
    std::pair <std::string, std::vector<double> > token( Rcpp::as<string>(names(i)), vec );
    map.insert( token );
  }
  return map;
}


// [[Rcpp::export]]
List getParameters_from_files(Rcpp::String folder)
{
    ecomeristem::ModelParameters parameters;
    utils::ParametersReader reader;
    reader.loadParametersFromFiles(folder, parameters);

    std::map < std::string, double > * paramMap = parameters.getRawParameters();
    Rcpp::List result(paramMap->size());
    Rcpp::CharacterVector names;
    Rcpp::NumericVector values;
    for(auto const& it: *paramMap){
        std::string key = it.first;
        double value = it.second;
        names.push_back(key);
        values.push_back(value);
    }

    DataFrame df = DataFrame::create(Named("Name")=names,Named("Values")=values);
    return df;
}


// [[Rcpp::export]]
List getMeteo_from_files(Rcpp::String folder)
{
    ecomeristem::ModelParameters parameters;
    utils::ParametersReader reader;
    reader.loadParametersFromFiles(folder, parameters);

    std::vector < ecomeristem::Climate > * meteoValues = parameters.getMeteoValues();
    Rcpp::List result(meteoValues->size());
    CharacterVector names =  CharacterVector::create("Temperature", "Par", "Etp", "Irrigation", "P");
    NumericVector Temperature, Par, Etp, Irrigation, P;

    for(auto const& it: *meteoValues){
        Temperature.push_back(it.Temperature);
        Par.push_back(it.Par);
        Etp.push_back(it.Etp);
        Irrigation.push_back(it.Irrigation);
        P.push_back(it.P);
    }

    DataFrame df = DataFrame::create(
                Named("Temperature")=Temperature,
                Named("Par")=Par,
                Named("Etp")=Etp,
                Named("Irrigation")=Irrigation,
                Named("P")=P
            );
    return df;
}


struct Simulation {
  Simulation(){}
  GlobalParameters globalParameters;
  ecomeristem::ModelParameters parameters;
  double beginDate;
  double endDate;
  EcomeristemContext context;
  SimulatorFilter filter;
};


std::map <std::string, Simulation*> simulations;

// [[Rcpp::export]]
void init_simu(List dfParameters, List dfMeteo, List obs, Rcpp::String name) {
  Simulation * s = new Simulation();
  simulations[name] = s;
  CharacterVector names = dfParameters[0];
  NumericVector values = dfParameters[1];
  for (int i = 0; i < names.size(); ++i) {
    s->parameters.getRawParameters()->insert(
        std::pair<std::string,double>(
            Rcpp::as<std::string>(names(i)),
            values(i)
        )
    );
  }

  NumericVector Temperature = dfMeteo[0];
  NumericVector Par = dfMeteo[1];
  NumericVector Etp = dfMeteo[2];
  NumericVector Irrigation = dfMeteo[3];
  NumericVector P = dfMeteo[4];

  for (int i = 0; i < Temperature.size(); ++i) {
    ecomeristem::Climate c(Temperature(i), Par(i), Etp(i), Irrigation(i), P(i));
    s->parameters.meteoValues.push_back(c);
  }

  s->parameters.beginDate = s->parameters.get("BeginDate");
  /** RUN SIMU **/
  s->beginDate = s->parameters.get("BeginDate");
  s->endDate = s->parameters.get("EndDate");
  s->context.setBegin(s->beginDate);
  s->context.setEnd(s->endDate);

  observer::PlantView view;
  map <string, vector <double > > obsMap = mapFromDF(obs);
  s->filter.init(&view, obsMap, "day");
}

// [[Rcpp::export]]
List launch_simu(Rcpp::String name, CharacterVector names = CharacterVector(), NumericVector params = NumericVector()) {
  Simulation * s = simulations[name];
  if(names.size() > 0) {
    for (int i = 0; i < names.size(); ++i) {
      s->parameters.mParams[Rcpp::as<string>(names(i))] = params[i];
    }
  }
  EcomeristemSimulator simulator(new PlantModel(), s->globalParameters);
  simulator.init(s->beginDate, s->parameters);
  map<string,vector<double>> res = simulator.runOptim(s->context, s->filter);
  return mapOfVectorToDF(res);
}

// [[Rcpp::export]]
List get_clean_obs(Rcpp::String vObsPath) {
  utils::ParametersReader reader;
  observer::PlantView view;
  std::map <std::string, std::vector<double> > vobsMap = reader.loadCleanObsFromFile(vObsPath, view);
  return mapOfVectorToDF(vobsMap);
}






// [[Rcpp::export]]
List rcpp_run_from_dataframe(List dfParameters, List dfMeteo)
{
  /** INIT PARAMS **/
  GlobalParameters globalParameters;
  ecomeristem::ModelParameters parameters;

  CharacterVector names = dfParameters[0];
  NumericVector values = dfParameters[1];
  for (int i = 0; i < names.size(); ++i) {
    parameters.getRawParameters()->insert(
        std::pair<std::string,double>(
            Rcpp::as<std::string>(names(i)),
            values(i)
        )
    );
  }

  NumericVector Temperature = dfMeteo[0];
  NumericVector Par = dfMeteo[1];
  NumericVector Etp = dfMeteo[2];
  NumericVector Irrigation = dfMeteo[3];
  NumericVector P = dfMeteo[4];

  for (int i = 0; i < Temperature.size(); ++i) {
    ecomeristem::Climate c(Temperature(i), Par(i), Etp(i), Irrigation(i), P(i));
    parameters.meteoValues.push_back(c);
  }

  parameters.beginDate = parameters.get("BeginDate");
  /** RUN SIMU **/
  double begin = parameters.get("BeginDate");
  double end = parameters.get("EndDate");

  EcomeristemSimulator * simulator = new EcomeristemSimulator(new PlantModel(), globalParameters);
  observer::PlantView *view = new observer::PlantView();
  simulator->attachView("plant", view);
  simulator->init(begin, parameters);
  EcomeristemContext context(begin, end);
  simulator->run(context);
  ResultParser parser;
  std::map <std::string, std::vector<double> > resultMap = parser.resultsToMap(simulator);
  delete simulator;
  return mapOfVectorToDF(resultMap);
}


// [[Rcpp::export]]
List rcpp_get_vObs(Rcpp::String vObsPath) {
  utils::ParametersReader reader;
  std::map <std::string, std::vector<double> > vobsMap = reader.loadVObsFromFile(vObsPath);
  return mapOfVectorToDF(vobsMap);
}

// [[Rcpp::export]]
List rcpp_reduceVobs(List vObs, List results) {
  std::map <std::string, std::vector<double> > vObsMap;
  std::map <std::string, std::vector<double> > resultMap;
  vObsMap = mapFromDF(vObs);
  resultMap = mapFromDF(results);
  ResultParser parser;
  auto ret = parser.filterVObs(vObsMap,resultMap, false);
  return mapOfVectorToDF(ret);
}


// [[Rcpp::export]]
List rcpp_reduceResults(List results, List vobs) {
  std::map <std::string, std::vector<double> > vObsMap;
  std::map <std::string, std::vector<double> > resultMap;
  vObsMap = mapFromDF(vobs);
  resultMap = mapFromDF(results);
  ResultParser parser;
  auto ret = parser.reduceResults(resultMap,vObsMap);
  return mapOfVectorToDF(ret);
}

