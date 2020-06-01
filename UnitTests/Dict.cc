// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdOdIdOdOdIKinKaldIUnitTestsdIDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../../KinKal/UnitTests/KKHitInfo.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *KinKalcLcLKKHitInfo_Dictionary();
   static void KinKalcLcLKKHitInfo_TClassManip(TClass*);
   static void *new_KinKalcLcLKKHitInfo(void *p = 0);
   static void *newArray_KinKalcLcLKKHitInfo(Long_t size, void *p);
   static void delete_KinKalcLcLKKHitInfo(void *p);
   static void deleteArray_KinKalcLcLKKHitInfo(void *p);
   static void destruct_KinKalcLcLKKHitInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::KinKal::KKHitInfo*)
   {
      ::KinKal::KKHitInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::KinKal::KKHitInfo));
      static ::ROOT::TGenericClassInfo 
         instance("KinKal::KKHitInfo", "../../KinKal/UnitTests/KKHitInfo.hh", 5,
                  typeid(::KinKal::KKHitInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &KinKalcLcLKKHitInfo_Dictionary, isa_proxy, 4,
                  sizeof(::KinKal::KKHitInfo) );
      instance.SetNew(&new_KinKalcLcLKKHitInfo);
      instance.SetNewArray(&newArray_KinKalcLcLKKHitInfo);
      instance.SetDelete(&delete_KinKalcLcLKKHitInfo);
      instance.SetDeleteArray(&deleteArray_KinKalcLcLKKHitInfo);
      instance.SetDestructor(&destruct_KinKalcLcLKKHitInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::KinKal::KKHitInfo*)
   {
      return GenerateInitInstanceLocal((::KinKal::KKHitInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::KinKal::KKHitInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *KinKalcLcLKKHitInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::KinKal::KKHitInfo*)0x0)->GetClass();
      KinKalcLcLKKHitInfo_TClassManip(theClass);
   return theClass;
   }

   static void KinKalcLcLKKHitInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_KinKalcLcLKKHitInfo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::KinKal::KKHitInfo : new ::KinKal::KKHitInfo;
   }
   static void *newArray_KinKalcLcLKKHitInfo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::KinKal::KKHitInfo[nElements] : new ::KinKal::KKHitInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_KinKalcLcLKKHitInfo(void *p) {
      delete ((::KinKal::KKHitInfo*)p);
   }
   static void deleteArray_KinKalcLcLKKHitInfo(void *p) {
      delete [] ((::KinKal::KKHitInfo*)p);
   }
   static void destruct_KinKalcLcLKKHitInfo(void *p) {
      typedef ::KinKal::KKHitInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::KinKal::KKHitInfo

namespace ROOT {
   static TClass *vectorlEKinKalcLcLKKHitInfogR_Dictionary();
   static void vectorlEKinKalcLcLKKHitInfogR_TClassManip(TClass*);
   static void *new_vectorlEKinKalcLcLKKHitInfogR(void *p = 0);
   static void *newArray_vectorlEKinKalcLcLKKHitInfogR(Long_t size, void *p);
   static void delete_vectorlEKinKalcLcLKKHitInfogR(void *p);
   static void deleteArray_vectorlEKinKalcLcLKKHitInfogR(void *p);
   static void destruct_vectorlEKinKalcLcLKKHitInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<KinKal::KKHitInfo>*)
   {
      vector<KinKal::KKHitInfo> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<KinKal::KKHitInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<KinKal::KKHitInfo>", -2, "vector", 469,
                  typeid(vector<KinKal::KKHitInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEKinKalcLcLKKHitInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<KinKal::KKHitInfo>) );
      instance.SetNew(&new_vectorlEKinKalcLcLKKHitInfogR);
      instance.SetNewArray(&newArray_vectorlEKinKalcLcLKKHitInfogR);
      instance.SetDelete(&delete_vectorlEKinKalcLcLKKHitInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlEKinKalcLcLKKHitInfogR);
      instance.SetDestructor(&destruct_vectorlEKinKalcLcLKKHitInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<KinKal::KKHitInfo> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<KinKal::KKHitInfo>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEKinKalcLcLKKHitInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<KinKal::KKHitInfo>*)0x0)->GetClass();
      vectorlEKinKalcLcLKKHitInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEKinKalcLcLKKHitInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEKinKalcLcLKKHitInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<KinKal::KKHitInfo> : new vector<KinKal::KKHitInfo>;
   }
   static void *newArray_vectorlEKinKalcLcLKKHitInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<KinKal::KKHitInfo>[nElements] : new vector<KinKal::KKHitInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEKinKalcLcLKKHitInfogR(void *p) {
      delete ((vector<KinKal::KKHitInfo>*)p);
   }
   static void deleteArray_vectorlEKinKalcLcLKKHitInfogR(void *p) {
      delete [] ((vector<KinKal::KKHitInfo>*)p);
   }
   static void destruct_vectorlEKinKalcLcLKKHitInfogR(void *p) {
      typedef vector<KinKal::KKHitInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<KinKal::KKHitInfo>

namespace {
  void TriggerDictionaryInitialization_Dict_Impl() {
    static const char* headers[] = {
"../../KinKal/UnitTests/KKHitInfo.hh",
0
    };
    static const char* includePaths[] = {
"/usr/local/Cellar/root/6.20.04_1/include/root",
"/Users/soleti/KFTrack/build_prof/UnitTests/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace KinKal{struct __attribute__((annotate("$clingAutoload$../../KinKal/UnitTests/KKHitInfo.hh")))  KKHitInfo;}
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "../../KinKal/UnitTests/KKHitInfo.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"KinKal::KKHitInfo", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dict() {
  TriggerDictionaryInitialization_Dict_Impl();
}
