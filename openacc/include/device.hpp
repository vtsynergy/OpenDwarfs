#ifndef __CTSAR_DEVICE_HPP
#define __CTSAR_DEVICE_HPP

struct ctsar;

namespace CTSAR{

class Device{
    private:
        //Scheduling master
        ctsar * master_;

        //Assignment information
        uint64_t assigned_;
        uint64_t start_;
        uint64_t end_;

        //Pass information, timing and iteration count
        struct timeval pass_start_time_;
        struct timeval pass_end_time_;
        double pass_time_; //time taken in this pass
        uint64_t pass_iterations_; //iterations completed by this device this pass

    public:
    Device() = delete;
    Device(ctsar *c)
         : master_(c){
    }
    virtual void assign(uint64_t start, uint64_t end){
        start_ = std::min(start, end);
        end_ = std::max(start,end);
        pass_iterations_ += start_ - end_;
    }
    virtual void start(){//begin pass, take timing
        gettimeofday(&pass_start_time_, NULL);
    }
    virtual void end(){//end pass, take timing
        gettimeofday(&pass_end_time_, NULL);
    }
    virtual void * allocate(uint64_t size){//allocate device memory
        throw "unimplemented";
    }
    virtual void copy_to(void * dest, void * src, uint64_t size){//copy into device space
        throw "unimplemented";
    }
    virtual void copy_from(void * dest, void * src, uint64_t size){//copy from device space
        throw "unimplemented";
    }

    virtual ~Device(){};
};
class DevNULL : public Device {//filler
    public:
    DevNULL(ctsar *c) : Device(c){
        std::cerr << "creating NULL device " << std::endl;
    }
    virtual void assign(uint64_t start,
                        uint64_t end){
        throw "unimplemented";
    }
};
class DevCPU : public Device {
    public:
    DevCPU(ctsar *c) : Device(c){
        std::cerr << "creating CPU " << std::endl;
    }
    virtual void assign(uint64_t start,
                        uint64_t end){
        throw "unimplemented";
    }
};
class DevGPU : public Device {
    public:
    DevGPU(ctsar * c) : Device(c){
        std::cerr << "creating GPU " << std::endl;
    }
    virtual void assign(uint64_t start,
                        uint64_t end){
        throw "unimplemented";
    }
};
class DevGPU_CUDA : public DevGPU {
    public:
    DevGPU_CUDA(ctsar * c) : DevGPU(c){
        std::cerr << "creating CUDA GPU " << std::endl;
    }
};
//aggregate device type, holds and schedules other devices
class DevGroup : public Device {
    public:
        DevGroup() = delete;
        DevGroup(ctsar * c, std::vector<std::shared_ptr<Device>> devices)
            : Device(c),
            devices_(devices) {
                group_barrier_ = std::make_shared<boost::barrier>(devices_.size());//here to ensure defined behavior of destroy
            }
        virtual ~DevGroup(){}

        //add another device, not recommended
        void add(std::shared_ptr<Device> dev){
            devices_.push_back(dev);
            group_barrier_ = std::make_shared<boost::barrier>(devices_.size());
        }
        virtual void assign(uint64_t start,
                            uint64_t end){
            throw "unimplemented";
        }
        virtual void start(){//begin pass, take timing
            throw "unimplemented";
        }
        virtual void end(){//end pass, take timing
            throw "unimplemented";
        }
    private:
        std::vector<std::shared_ptr<Device>> devices_;
        std::shared_ptr<boost::barrier> group_barrier_;
};
}

#endif // __CTSAR_DEVICE_HPP
