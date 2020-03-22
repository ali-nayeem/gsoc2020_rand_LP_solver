
#include "Eigen/Eigen"
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include "random.hpp"
#include "random/uniform_int.hpp"
#include "random/normal_distribution.hpp"
#include "random/uniform_real_distribution.hpp"
#include "volume.h"
#include <hpolytope.h>
#include <fstream>
#include <unordered_map>
#include <iostream>
#include <numeric>
//#define VOLESTI_DEBUG

using namespace std;

template <class Point>
class lp_problem
{
    unordered_map<string, unsigned> decisionVar;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    typedef boost::mt19937 RNGType;
    RNGType rng; //seed
public:
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    string direction;
    HPolytope<Point> polytope;
    VT objectiveFunction;
    Point solution;
    NT solutionVal;
    int dim;
    int row;
    Point obj;

    /* lp_solve LP file format
    TODO: Consider comment, blank line, coefficient without sign
    */
    lp_problem(string pathToProblem)
    {
        seed = std::chrono::system_clock::now().time_since_epoch().count();
        rng =  RNGType(seed); //seed
        std::ifstream file(pathToProblem, std::ifstream::in);

        std::string line;
        std::getline(file, line, '\n');
        if (line.at(line.size() - 1) != ';')
        {
            throw "line not ending with ;";
        }
        line.resize(line.size() - 1); //eliminate ; at the end of each line
        std::istringstream objStream(line);
        std::string aWord;
        objStream >> aWord;
        direction = aWord; //min:
        dim = 0;
        while (objStream >> aWord) //scanning the objective function
        {
            objectiveFunction.conservativeResize(dim + 1);
            try
            {
                double coef = std::stod(aWord);
                objectiveFunction(dim) = coef;
                objStream >> aWord; //get the name of the variable
                decisionVar[aWord] = dim;
            }
            catch (const std::invalid_argument &ia)
            {
                /* ... */
                objectiveFunction(dim) = (aWord.at(0) == '+') ? 1.0 : -1.0;
                aWord = aWord.substr(1); //eliminate sign
                decisionVar[aWord] = dim;
            }
            dim++;
            if (direction.find("min") == std::string::npos) //not min => max
            {
                objectiveFunction(dim - 1) = -1 * objectiveFunction(dim - 1);
            }
        }
        //read constraints
        std::list<std::vector<NT>> constraints;
        VT b;
        row = 0;
        while (!std::getline(file, line, '\n').eof())
        {
            line.resize(line.size() - 1); //eliminate ; at the end of each line
            std::istringstream constraintStream(line);
            constraintStream >> aWord;
            std::vector<NT> aConstraint(dim);
            while (constraintStream >> aWord) // for each constraint
            {
                if (aWord.find("<") != std::string::npos || aWord.find(">") != std::string::npos)
                {
                    string relation = aWord;
                    constraintStream >> aWord; //rhs
                    b.conservativeResize(row + 1);
                    b(row) = std::stod(aWord);
                    if (relation.find(">") != std::string::npos) //converting a GT constraint into a LT
                    {
                        b(row) = -1 * b(row);
                        for (int i = 0; i < aConstraint.size(); i++)
                        {
                            aConstraint[i] = -1 * aConstraint[i];
                        }
                    }

                    row++;
                    break;
                }

                try
                {
                    double coef = std::stod(aWord);
                    // objectiveFunction(dim) = coef;
                    constraintStream >> aWord; //get the name of the variable
                    // decisionVar[aWord] = dim++;
                    if (decisionVar.find(aWord) == decisionVar.end())
                    {
                        aConstraint.resize(dim + 1, 0);
                        objectiveFunction.resize(dim + 1);
                        objectiveFunction(dim) = 0;
                        decisionVar[aWord] = dim++;
                    }
                    aConstraint[decisionVar[aWord]] = coef;
                }
                catch (const std::invalid_argument &ia)
                {
                    double coef = (aWord.at(0) == '+') ? 1.0 : -1.0;
                    aWord = aWord.substr(1); //eliminate sign
                    //decisionVar[aWord] = dim++;
                    if (decisionVar.find(aWord) == decisionVar.end())
                    {
                        aConstraint.resize(dim + 1, 0);
                        objectiveFunction.resize(dim + 1);
                        objectiveFunction(dim) = 0;
                        decisionVar[aWord] = dim++;
                    }
                    aConstraint[decisionVar[aWord]] = coef;
                }
            }
            constraints.push_back(aConstraint);
        }

        MT A;
        A.resize(row, dim);

        int i = 0;
        for (auto constraint : constraints)
        {
            for (int j = 0; j < dim; j++)
            {
                if (j < constraint.size())
                {
                    A(i, j) = constraint[j];
                }
                else
                {
                    A(i, j) = 0;
                }
            }
            i++;
        }
        this->polytope.init(dim, A, b);
        //this->polytope.print();
        this->obj = Point(this->dim);
        for (int i = 0; i < this->dim; i++)
        {
            obj.set_coord(i, objectiveFunction(i));
        }
    }


    Point learnInteriorPoint(MT A, VT b, double learningRate = 0.1, int numOfEpoch = 1000)
    {
        Point weights(this->dim, vector<double>(this->dim, 1)); 
        // init the walks
        //RNGType &rng = var.rng;
        boost::random::uniform_real_distribution<> urdist(-1, 1);
        for (size_t i = 0; i < this->dim; i++)
        {
            weights.set_coord(i, urdist(rng));
        }
        double minError = std::numeric_limits<double>::max();
        Point bestWeight(this->dim);
        for (int epoch = 0; epoch < numOfEpoch; epoch++) // for epoch in range(n_epoch):
        {

            double sum_error = 0.0;
            for (size_t r = 0; r < row; r++) //for row in train:
            {
                Point row(this->dim);
                for (size_t i = 0; i < this->dim; i++)
                {
                    row.set_coord(i, A(r, i));
                }
                double prediction = weights.dot(row);
                double error = b(r) - prediction;
                error = (error >= 0) ? 0 : error;
                sum_error += error * -1;
                for (size_t i = 0; i < this->dim; i++)
                {
                    weights.set_coord(i, weights[i] + learningRate * error * row[i]);
                }
            }

            if (sum_error < minError)
            {
                minError = sum_error;
                bestWeight = Point(dim, weights.get_coeffs());
            }
            printf(">epoch=%d, lrate=%.3f, error=%f\n", epoch, learningRate, sum_error);
            if (sum_error < 0.00001)
            {
                polytope.is_in(weights);
                break;
            }
        }
        return bestWeight;
    }
    void rand_point_generator(HPolytope<Point> &P, Point &p, const unsigned int rnum, const unsigned int walk_len, std::vector<Point> &randPoints)
    {
        randPoints.clear();
        // init the walks
        boost::random::uniform_real_distribution<> urdist(0, 1);
        boost::random::uniform_int_distribution<> uidist(0, dim - 1);
        std::vector<NT> lamdas(P.num_of_hyperplanes(), NT(0));
        unsigned int rand_coord, rand_coord_prev;
        NT kapa;
        Point p_prev = p;

        //Compute the first point for the CDHR
        rand_coord = uidist(rng);
        kapa = urdist(rng);

        std::pair<NT, NT> bpair = P.line_intersect_coord(p, rand_coord, lamdas);
        p.set_coord(rand_coord, p[rand_coord] + bpair.first + kapa * (bpair.second - bpair.first));

        for (unsigned int i = 1; i <= rnum; ++i)
        {
            for (unsigned int j = 0; j < walk_len; ++j)
            {
                rand_coord_prev = rand_coord;
                rand_coord = uidist(rng);
                kapa = urdist(rng);
                hit_and_run_coord_update(p, p_prev, P, rand_coord, rand_coord_prev, kapa, lamdas);
            }
            if (P.is_in(p) != 0)
            {
                randPoints.push_back(p);
            }
        }
    }
    void randCuttingPlane(int maxIter=1000)
    {
        //polytope.print();
        MT A = polytope.get_mat();
        VT b = polytope.get_vec();
        A.conservativeResize(this->row + 1, dim);
        b.conservativeResize(this->row + 1);
        for (int i = 0; i < dim; i++)
        {
            A(row, i) = objectiveFunction(i);
        }

        std::vector<Point> randPoints;
        Point minPoint(this->dim);
        NT minObj;
        for (int iter = 0; iter < maxIter; iter++)
        {
            std::pair<Point, NT> InnerBall = polytope.ComputeInnerBall();
 
            Point start = InnerBall.first; //learnInteriorPoint(polytope.get_mat(), polytope.get_vec());
            rand_point_generator(polytope, start, 10, 3, randPoints);
            int c = 0;
            for (auto pt : randPoints)
            {
                if (polytope.is_in(pt) == 0)
                {
                    c++;
                }
            }
            //cout << c << " Not in poly\n";
            
            minPoint = randPoints[0];
            minObj = minPoint.dot(this->obj);
            vector<NT> objValues;
            objValues.push_back(minObj);
            for (size_t i = 1; i < randPoints.size(); i++)
            {
                NT val = randPoints[i].dot(this->obj);
                objValues.push_back(val);
                if (val < minObj)
                {
                    minObj = val;
                    minPoint = randPoints[i];
                }
            }

            double sum = std::accumulate(objValues.begin(), objValues.end(), 0.0);
            double mean = sum / objValues.size();
            std::vector<double> diff(objValues.size());
            std::transform(objValues.begin(), objValues.end(), diff.begin(), [mean](double x) { return x - mean; });
            double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / objValues.size());
            if (stdev < 0.00001)
            {
                break;
            }
            b(row) = minObj;
            polytope.init(dim, A, b);
            // polytope.set_mat(A);
            // polytope.set_vec(b);
            cout << "iteration: " << iter << ", objective: " << minObj<<endl;
        }
        solution = minPoint;
        solutionVal = minObj;
    }
    void printSolution()
    {
        std::cout << "Best score is: " << solutionVal << std::endl << "coords: ";
        solution.print();
    }
};

int main()
{
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    lp_problem<Point> prob("test.lp"); //sc50a afiro kb2
    prob.randCuttingPlane();
    prob.printSolution();
    //prob.train_weights(prob.polytope.get_mat(), prob.polytope.get_vec());
    return 0;
}
//#endif //VOLESTI_LP_PROBLEM_H