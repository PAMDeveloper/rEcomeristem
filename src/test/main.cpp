// #include <ctime>
#define UNSAFE_RUN
#include <defines.hpp>
#include <utils/ParametersReader.hpp>
#include <plant/PlantModel.hpp>
#include <observer/PlantView.hpp>

#include <ctime>
using namespace std;

void display(map<string, vector<double>> map) {
	for (auto token : map) {
		//        if(token.first != "LIG" && token.first != "lig")
		//            continue;
		cout << "**************" << token.first << "**************" << "\n";
		string vals = "";
		for (double val : token.second) {
			vals += to_string(val) + ",";
		}
		cout << vals << "\n";
	}
}
#ifdef UNSAFE_RUN
//void traceModel(AbstractSimpleModel * m, double step, QString parent = "", double idx = -1) {
//    QString name = (parent.isEmpty() ? "" : "/")
//                    + (idx > -1 ? "[" + QString::number(idx) + "]" : "")
//                    + QString::fromStdString(m->name());

//    for (int i = 0; i < m->i_names.size(); ++i) {
//        qDebug() << fixed <<  step << " - " << name << " : "
//                 << QString::fromStdString(m->i_names[i])
//                 << m->getIVal(i, step, true)
//                 << " > " << m->getIVal(i, step, false);
//    }
//    for(auto token: m->subModels) {
//        for (int idx = 0; idx < token.second.size(); ++idx) {
//            traceModel(token.second[idx], step, name, idx);
//        }
//    }
//}
#endif

#include <ctime>
struct Simulation {
public:
	Simulation() {}
	GlobalParameters globalParameters;
	ecomeristem::ModelParameters parameters;
	double beginDate;
	double endDate;
	EcomeristemContext context;
	SimulatorFilter filter;
};

int main(int argc, char *argv[]) {

	std::string dirName = "D:\\Samples\\_Estimation\\G1";
	ecomeristem::ModelParameters parameters;
	utils::ParametersReader reader;
	reader.loadParametersFromFiles(dirName, parameters);
	std::map <std::string, std::vector<double> > obsMap = reader.loadVObsFromFile(dirName + "\\vobs_moy.txt");
	Simulation * s = new Simulation();
	s->parameters = parameters;
	s->parameters.beginDate = s->parameters.get("BeginDate");
	s->beginDate = s->parameters.get("BeginDate");
	s->endDate = s->parameters.get("EndDate");
	s->context.setBegin(s->beginDate);
	s->context.setEnd(s->endDate);
	observer::PlantView view;
	s->filter.init(&view, obsMap, "day");
	EcomeristemSimulator simulator(new PlantModel(), s->globalParameters);
	simulator.init(s->beginDate, s->parameters);
	map<string, vector<double>> res = simulator.runOptim(s->context, s->filter);
	display(res);
	display(obsMap);
	return 1;

    /***TIMER***/

 //   std::string dirName = "D:\\Samples\\_Estimation\\G1";
 //   ecomeristem::ModelParameters parameters;
 //   utils::ParametersReader reader;
 //   reader.loadParametersFromFiles(dirName, parameters);
 //   // qDebug() << fixed << parameters.beginDate << parameters.get("EndDate");
 //   const clock_t begin_time = clock();
 //   GlobalParameters globalParameters;
 //   EcomeristemContext context(parameters.get("BeginDate"), parameters.get("EndDate"));
	//observer::PlantView *view = new observer::PlantView();
 //   for(int i = 0; i < 1000; i++) {
 //       PlantModel * m = new PlantModel;
 //       EcomeristemSimulator simulator(m, globalParameters);
 //       simulator.attachView("plant", view);
 //       simulator.init(parameters.get("BeginDate"), parameters);
 //       simulator.run(context);
 //   }
	//delete view;

    //std::cout << float( clock () - begin_time ) /  (CLOCKS_PER_SEC) << " ms";
    return 1;
}