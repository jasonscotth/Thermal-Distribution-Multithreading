#include<iostream>
#include<iomanip>
#include<sys/time.h>
#include<pthread.h>
#include<cmath>

using namespace std;

//thread argument struct definition
struct thread_struct {
    //unique thread id
    int threadid;
    //local start and ending positions
    int start;
    int end;
    //copies of edge markers
    int edge_h;
    int edge_w;
    //pointing to a single mutex/barrier declared in main
    pthread_mutex_t *mutex;
    pthread_barrier_t *barrier;
};

//global variables
float global_max = 0;
float epsilon = 0.1;
float*** cube = nullptr;
float*** cube_update = nullptr;
thread_struct* thread_info = nullptr;
int* partitions;
int* start_indices;
int* end_indices;
bool global_flag = true;
int iteration = 0;
int interval = 1;

//allocate memory, initialize faces and internal points for data structures
void create_structures(int depth, int height, int width, int edge_d, int edge_h, int edge_w,
int front, int back, int left, int right, int top, int bottom) {
    cout << "Building data structures in memory..." << endl;

    //depth of cube
    cube = new float**[depth];
    cube_update = new float**[depth];

    //rows
    for(int i=0; i<depth; i++) {
        cube[i] = new float*[height];
        cube_update[i] = new float*[height];
    }

    //columns
    for(int i=0; i<depth; i++) {
        for(int j=0; j<height; j++) {
            cube[i][j] = new float[width];
            cube_update[i][j] = new float[width];
        }
    }

    cout << "Memory allocated!" << endl;
    cout << "Initializing edges..." << endl;

    //create the front and back faces
    for(int j=0; j<height; j++) {
        for(int k=0; k<width; k++) {
            cube[0][j][k] = front;
            cube_update[0][j][k] = front;
            cube[edge_d][j][k] = back;
            cube_update[edge_d][j][k] = back;
        }
    }

    //create the left and right faces,
    //then the top and bottom faces that overwrite left and right
    for(int i=1; i<edge_d; i++) {
        for(int j=1; j<edge_h; j++) {
            cube[i][j][0] = left;
            cube[i][j][edge_w] = right;
            cube_update[i][j][0] = left;
            cube_update[i][j][edge_w] = right;
        }
        for(int k=0; k<width; k++) {
            cube[i][0][k] = top;
            cube[i][edge_h][k] = bottom;
            cube_update[i][0][k] = top;
            cube_update[i][edge_h][k] = bottom;
        }
    }

    cout << "Edges initialized!" << endl;
    cout << "Initializing internal points..." << endl;

    //calculate internal point value
    //front/back dimensions
    float x = height * width;
    //top and bottom dimensions
    float y = width * (depth - 2);
    //left and right dimensions
    float z = (depth - 2) * (height - 2);

    float step1 = ((front * x) + (back * x)) / x;
    float step2 = ((top * y) + (bottom * y)) / y;
    float step3 = ((left * z) + (right * z)) / z;
    float initialize = (step1 + step2 + step3) / 3;

    cout << "internal points value: " << initialize << endl;

    for(int i=1; i<edge_d; i++) {
        for(int j=1; j<edge_h; j++) {
            for(int k=1; k<edge_w; k++) {
                cube[i][j][k] = initialize;
            }
        }
    }

    cout << "Internal points initialized!" << endl;
}

//find split points in structures for multiprocessing
void partition(int num_threads, int depth, int edge_h, int edge_w, pthread_mutex_t *mutex, pthread_barrier_t *barrier) {
    partitions = new int[num_threads];
    start_indices = new int[num_threads];
    end_indices = new int[num_threads];
    //find equal division of depth of structure
    int division = (depth-2) / num_threads;
    //find remainder, if any, of division
    int modulus = (depth-2) % num_threads;
    int start_pos = 1;
    int end_pos = 1;
    for(int i=0; i<num_threads; i++) {
        partitions[i] = division;
        if(modulus != 0){
            partitions[i] += 1;
            modulus -= 1;
        }
    }
    for(int i=0; i<num_threads; i++) {
        start_indices[i] = start_pos;
        end_pos += partitions[i];
        end_indices[i] = end_pos;
        start_pos += partitions[i];
    }

    thread_info = new thread_struct[num_threads];
    for(int i=0; i<num_threads; i++) {
        thread_info[i].threadid = i;
        thread_info[i].start = start_indices[i];
        thread_info[i].end = end_indices[i];
        thread_info[i].edge_h = edge_h;
        thread_info[i].edge_w = edge_w;
        thread_info[i].mutex = mutex;
        thread_info[i].barrier = barrier;
        //print test
        //cout << myStruct[i].threadid << " " << myStruct[i].start << " " << myStruct[i].end << endl;
    }
    //memory efficiency
    delete[] partitions;
    partitions = nullptr;
    delete[] start_indices;
    start_indices = nullptr;
    delete[] end_indices;
    end_indices = nullptr;
};

void *thread_work(void *thread_info) {
    thread_struct *arg_struct = (struct thread_struct *) thread_info;
    int id = arg_struct->threadid;
    int start = arg_struct->start;
    int end = arg_struct->end;
    int edge_h = arg_struct->edge_h;
    int edge_w = arg_struct->edge_w;

    while(global_flag) {
        //cout << "Starting calculations..." << endl;
        float local_max = 0;
        float temp = 0;
        int ret = 0;

        //cout << "Iteration: " << iteration << " by thread: " << id << endl;

        if(iteration % 2 == 0) {
            //cout << iteration << " Updating cube_update" << endl;
            for(int i=start; i<end; i++) {
                for(int j=1; j<edge_h; j++) {
                    for(int k=1; k<edge_w; k++) {
                        cube_update[i][j][k] = (cube[i+1][j][k] + cube[i-1][j][k] + cube[i][j+1][k]
                        + cube[i][j-1][k] + cube[i][j][k+1] + cube[i][j][k-1]) / 6;
                        temp = abs(cube_update[i][j][k] - cube[i][j][k]);
                        //cout << temp << endl;
                        if(temp > local_max) {local_max = temp;}
                    }
                }
            }
            if(local_max > global_max) {
                pthread_mutex_lock(arg_struct->mutex);
                global_max = local_max;
                pthread_mutex_unlock(arg_struct->mutex);
            }
        }
        else {
            //cout << iteration << " Updating cube" << endl;
            for(int i=start; i<end; i++) {
                for(int j=1; j<edge_h; j++) {
                    for(int k=1; k<edge_w; k++) {
                        cube[i][j][k] = (cube_update[i+1][j][k] + cube_update[i-1][j][k] + cube_update[i][j+1][k]
                        + cube_update[i][j-1][k] + cube_update[i][j][k+1] + cube_update[i][j][k-1]) / 6;
                        temp = abs(cube[i][j][k] - cube_update[i][j][k]);
                        //cout << temp << endl;
                        if(temp > local_max) {local_max = temp;}
                    }
                }
            }
            if(local_max > global_max) {
                pthread_mutex_lock(arg_struct->mutex);
                global_max = local_max;
                pthread_mutex_unlock(arg_struct->mutex);
            }
        }
        ret = pthread_barrier_wait(arg_struct->barrier);
        ret = pthread_barrier_wait(arg_struct->barrier);
    }

    pthread_exit(NULL);
}

int main() {
    //define number of threads to be used
    int num_threads = 8;

    //dimension arguments
    int depth = 100;
    int height = 100;
    int width = 100;
    int edge_d = depth-1;
    int edge_h = height-1;
    int edge_w = width-1;

    //face temperatures
    float front = 600;
    float back = 800;
    float left = 400;
    float right = 500;
    float top = 100;
    float bottom = 100;

    //structs to hold runtime markers
    struct timeval start, end;
    gettimeofday(&start,NULL);

    //output run info
    cout << "Dimensions for your cube temperature simulation: " << endl;
    cout << "Depth: " << depth << "\nHeight: " << height << "\nWidth: " << width << endl;

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, num_threads);
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex,NULL);

    //allocate memory, initialize faces and internal points for data structures
    create_structures(depth,height,width,edge_d,edge_h,edge_w,front,back,left,right,top,bottom);

    //find split points in structures for multiprocessing
    partition(num_threads, depth, edge_h, edge_w, &mutex, &barrier);

    //thread creation and struct argument passing
    pthread_t threads[num_threads];
    //create 7 child threads, with main still participating in work
    for(int i=1; i<num_threads; i++) {
        //cout << "Creating thread " << i << endl;
        //return code from creation of POSIX thread given threadid, null, function to perform, and argument
        int rc = pthread_create(&threads[i], NULL, thread_work, (void *)&thread_info[i]);
        if (rc) {
            cout << "ERROR; return code from pthread_create() is" << rc << endl;
            exit(-1);
        }
    }

    while(global_flag) {
        //cout << "Starting calculations..." << endl;
        float local_max = 0;
        float temp = 0;
        int ret = 0;

        //cout << "Iteration: " << iteration << " by thread: 0" << endl;

        if(iteration % 2 == 0) {
            //cout << iteration << " Updating cube_update" << endl;
            for(int i=thread_info[0].start; i<thread_info[0].end; i++) {
                for(int j=1; j<edge_h; j++) {
                    for(int k=1; k<edge_w; k++) {
                        cube_update[i][j][k] = (cube[i+1][j][k] + cube[i-1][j][k] + cube[i][j+1][k]
                        + cube[i][j-1][k] + cube[i][j][k+1] + cube[i][j][k-1]) / 6;
                        temp = abs(cube_update[i][j][k] - cube[i][j][k]);
                        //cout << temp << endl;
                        if(temp > local_max) {local_max = temp;}
                    }
                }
            }
            if(local_max > global_max) {
                pthread_mutex_lock(&mutex);
                global_max = local_max;
                pthread_mutex_unlock(&mutex);
            }
        }
        else {
            //cout << iteration << " Updating cube" << endl;
            for(int i=thread_info[0].start; i<thread_info[0].end; i++) {
                for(int j=1; j<edge_h; j++) {
                    for(int k=1; k<edge_w; k++) {
                        cube[i][j][k] = (cube_update[i+1][j][k] + cube_update[i-1][j][k] + cube_update[i][j+1][k]
                        + cube_update[i][j-1][k] + cube_update[i][j][k+1] + cube_update[i][j][k-1]) / 6;
                        temp = abs(cube[i][j][k] - cube_update[i][j][k]);
                        //cout << temp << endl;
                        if(temp > local_max) {local_max = temp;}
                    }
                }
            }
            if(local_max > global_max) {
                pthread_mutex_lock(&mutex);
                global_max = local_max;
                pthread_mutex_unlock(&mutex);
            }
        }

        //first barrier, where all threads found their local max and global max has been determined
        ret = pthread_barrier_wait(&barrier);

        //cout << "The maximum difference found was: " << maximum_diff << endl;

        if(global_max < epsilon) {global_flag = false;}

        //print output
        //if flag set to false, print final interation
        if(!global_flag) {
            //print current iteration and end loop
            cout << setw(7) << iteration << setw(12) << fixed << setprecision(6) << global_max << endl;
        }
        //if flag set to true, print at 2^n intervals and continue
        else {
            if(iteration == interval) {
                cout << setw(7) << interval << setw(12) << fixed << setprecision(6) << global_max << endl;
                interval = interval << 1;
            }
        }
        iteration++;

        global_max = 0;

        //allow other threads to continue
        ret = pthread_barrier_wait(&barrier);
    }

    cout << "Freeing memory..." << endl;

    //depth/sheet iterator
    for(int i=0; i<depth; i++) {
        //row iterator
        for(int j=0; j<height; j++) {
            delete[] cube[i][j];
            delete[] cube_update[i][j];
        }
    }

    for(int i=0; i<depth; i++) {
        delete[] cube[i];
        delete[] cube_update[i];
    }

    delete[] cube;
    delete[] cube_update;

    cout << "Memory freed!" << endl;

    pthread_mutex_destroy(&mutex);
    pthread_barrier_destroy(&barrier);

    gettimeofday(&end,NULL);

    //calculate time difference
    double end_sec = end.tv_sec;
    double start_sec = start.tv_sec;
    double end_usec = end.tv_usec;
    double start_usec = start.tv_usec;
    double seconds_taken = (end_sec - start_sec) + ((end_usec - start_usec) / 1000000);
    int minutes_taken = seconds_taken / 60;
    int seconds = seconds_taken;
    seconds = seconds % 60;
    cout << minutes_taken << " Minutes " << seconds << "Seconds taken" << endl;

    return 0;
}