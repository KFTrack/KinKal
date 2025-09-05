// Ed Callaghan
// Memoization context for reinstantiating objects with graph relations
// August 2025

#ifndef KinKal_CloneContext_hh
#define KinKal_CloneContext_hh

#include <unordered_map>
#include <memory>

class CloneContext{
  public:
    // compiler defaults propagate to underlying stl-map
    CloneContext() = default;
   ~CloneContext() = default;
    // disallow copy-constructor so as to prevent bugs stemming from
    // failure to propagate instantiations across references to context
    CloneContext(CloneContext const&) = delete;

    // standard interface:
    // - if first call for domain ptr, instantiate a clone
    // - return ptr to clone
    template<typename T>
    std::shared_ptr<T> get(const std::shared_ptr<T>&);

    void clear();

  protected:
    // map of raw address to stl-pointers
    std::unordered_map <void*, std::shared_ptr<void> > map_;

  private:
    /**/
};

template<typename T>
std::shared_ptr<T> CloneContext::get(const std::shared_ptr<T>& ptr){
  void* address = static_cast<void*>(ptr.get());
  if (this->map_.count(address) < 1){
    auto cloned = ptr->clone(*this);
    this->map_[address] = static_pointer_cast<void>(cloned);
  }
  auto mapped = this->map_[address];
  auto rv = static_pointer_cast<T>(mapped);
  return rv;
}

#endif
