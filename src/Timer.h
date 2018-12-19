class Timer
{
public:

  // Constructor and distructor

  Timer();
  ~Timer();	

  void start();
  void reset();
  double stop();
  void check();

  double CPUTime();
  
private:
  clock_t 	_startCPUTime;	// start CPU time
  double	_startOMPTime;	// start CPU time (for OpenMP)
  double	_sumCPUtime;	// elapsed CPU time in seconds
  bool 		_isTimerOn;	// true if timer is running, false elswhere
  
};

