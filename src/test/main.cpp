// #include <ctime>
#define UNSAFE_RUN
#include <defines.hpp>
#include <utils/ParametersReader.hpp>
#include <plant/PlantModel.hpp>
#include <observer/PlantView.hpp>

#include <ctime>
using namespace std;

int main(int argc, char *argv[]) {

    /***TIMER***/

    std::string dirName = "D:\\Samples\\_Estimation\\G1";
    ecomeristem::ModelParameters parameters;
    utils::ParametersReader reader;
    reader.loadParametersFromFiles(dirName, parameters);
    // qDebug() << fixed << parameters.beginDate << parameters.get("EndDate");
    const clock_t begin_time = clock();
    GlobalParameters globalParameters;
    EcomeristemContext context(parameters.get("BeginDate"), parameters.get("EndDate"));
    for(int i = 0; i < 100; i++) {
        PlantModel * m = new PlantModel;
        EcomeristemSimulator simulator(m, globalParameters);
        observer::PlantView *view = new observer::PlantView();
        simulator.attachView("plant", view);
        simulator.init(parameters.get("BeginDate"), parameters);
        simulator.run(context);
        delete view;
    }
    std::cout << float( clock () - begin_time ) /  (CLOCKS_PER_SEC * 100);
    return 1;
}