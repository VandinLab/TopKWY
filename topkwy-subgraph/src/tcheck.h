extern Frequency MINFREQ;
extern int N_SMALL, N_TOTAL;
extern double ALPHA, THRESHOLD, COUNT, DELTA;
extern double wy_processing_time;
extern vector<int> CLASS_VEC;
extern vector<Tid> OCC_VEC;

void readClass(char *class_file);
void outputSubgraph(int count, Frequency frequency, double p_value);
void checkTestable(vector<LegOccurrence>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info);
void checkTestableCl(vector<CloseLegOccurrence>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info);
void checkCondition();
template<typename T> void checkConditionWY(vector<T>& elements, Frequency frequency, bounds_info& parent_info, bounds_info& current_info);
template<typename T> void computePvalue(vector<T>& elements, Frequency frequency);
