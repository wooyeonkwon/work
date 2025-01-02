// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME classes
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "SelectedMuon.h"
#include "ZBosonInfo.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *SelectedMuon_Dictionary();
   static void SelectedMuon_TClassManip(TClass*);
   static void *new_SelectedMuon(void *p = nullptr);
   static void *newArray_SelectedMuon(Long_t size, void *p);
   static void delete_SelectedMuon(void *p);
   static void deleteArray_SelectedMuon(void *p);
   static void destruct_SelectedMuon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SelectedMuon*)
   {
      ::SelectedMuon *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SelectedMuon));
      static ::ROOT::TGenericClassInfo 
         instance("SelectedMuon", "SelectedMuon.h", 10,
                  typeid(::SelectedMuon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SelectedMuon_Dictionary, isa_proxy, 4,
                  sizeof(::SelectedMuon) );
      instance.SetNew(&new_SelectedMuon);
      instance.SetNewArray(&newArray_SelectedMuon);
      instance.SetDelete(&delete_SelectedMuon);
      instance.SetDeleteArray(&deleteArray_SelectedMuon);
      instance.SetDestructor(&destruct_SelectedMuon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SelectedMuon*)
   {
      return GenerateInitInstanceLocal(static_cast<::SelectedMuon*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::SelectedMuon*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SelectedMuon_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::SelectedMuon*>(nullptr))->GetClass();
      SelectedMuon_TClassManip(theClass);
   return theClass;
   }

   static void SelectedMuon_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ZBosonInfo_Dictionary();
   static void ZBosonInfo_TClassManip(TClass*);
   static void *new_ZBosonInfo(void *p = nullptr);
   static void *newArray_ZBosonInfo(Long_t size, void *p);
   static void delete_ZBosonInfo(void *p);
   static void deleteArray_ZBosonInfo(void *p);
   static void destruct_ZBosonInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ZBosonInfo*)
   {
      ::ZBosonInfo *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ZBosonInfo));
      static ::ROOT::TGenericClassInfo 
         instance("ZBosonInfo", "ZBosonInfo.h", 5,
                  typeid(::ZBosonInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ZBosonInfo_Dictionary, isa_proxy, 4,
                  sizeof(::ZBosonInfo) );
      instance.SetNew(&new_ZBosonInfo);
      instance.SetNewArray(&newArray_ZBosonInfo);
      instance.SetDelete(&delete_ZBosonInfo);
      instance.SetDeleteArray(&deleteArray_ZBosonInfo);
      instance.SetDestructor(&destruct_ZBosonInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ZBosonInfo*)
   {
      return GenerateInitInstanceLocal(static_cast<::ZBosonInfo*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ZBosonInfo*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ZBosonInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::ZBosonInfo*>(nullptr))->GetClass();
      ZBosonInfo_TClassManip(theClass);
   return theClass;
   }

   static void ZBosonInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_SelectedMuon(void *p) {
      return  p ? new(p) ::SelectedMuon : new ::SelectedMuon;
   }
   static void *newArray_SelectedMuon(Long_t nElements, void *p) {
      return p ? new(p) ::SelectedMuon[nElements] : new ::SelectedMuon[nElements];
   }
   // Wrapper around operator delete
   static void delete_SelectedMuon(void *p) {
      delete (static_cast<::SelectedMuon*>(p));
   }
   static void deleteArray_SelectedMuon(void *p) {
      delete [] (static_cast<::SelectedMuon*>(p));
   }
   static void destruct_SelectedMuon(void *p) {
      typedef ::SelectedMuon current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::SelectedMuon

namespace ROOT {
   // Wrappers around operator new
   static void *new_ZBosonInfo(void *p) {
      return  p ? new(p) ::ZBosonInfo : new ::ZBosonInfo;
   }
   static void *newArray_ZBosonInfo(Long_t nElements, void *p) {
      return p ? new(p) ::ZBosonInfo[nElements] : new ::ZBosonInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_ZBosonInfo(void *p) {
      delete (static_cast<::ZBosonInfo*>(p));
   }
   static void deleteArray_ZBosonInfo(void *p) {
      delete [] (static_cast<::ZBosonInfo*>(p));
   }
   static void destruct_ZBosonInfo(void *p) {
      typedef ::ZBosonInfo current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ZBosonInfo

namespace ROOT {
   static TClass *vectorlEZBosonInfogR_Dictionary();
   static void vectorlEZBosonInfogR_TClassManip(TClass*);
   static void *new_vectorlEZBosonInfogR(void *p = nullptr);
   static void *newArray_vectorlEZBosonInfogR(Long_t size, void *p);
   static void delete_vectorlEZBosonInfogR(void *p);
   static void deleteArray_vectorlEZBosonInfogR(void *p);
   static void destruct_vectorlEZBosonInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ZBosonInfo>*)
   {
      vector<ZBosonInfo> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ZBosonInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ZBosonInfo>", -2, "vector", 423,
                  typeid(vector<ZBosonInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEZBosonInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ZBosonInfo>) );
      instance.SetNew(&new_vectorlEZBosonInfogR);
      instance.SetNewArray(&newArray_vectorlEZBosonInfogR);
      instance.SetDelete(&delete_vectorlEZBosonInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlEZBosonInfogR);
      instance.SetDestructor(&destruct_vectorlEZBosonInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ZBosonInfo> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ZBosonInfo>","std::vector<ZBosonInfo, std::allocator<ZBosonInfo> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ZBosonInfo>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEZBosonInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ZBosonInfo>*>(nullptr))->GetClass();
      vectorlEZBosonInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEZBosonInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEZBosonInfogR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ZBosonInfo> : new vector<ZBosonInfo>;
   }
   static void *newArray_vectorlEZBosonInfogR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ZBosonInfo>[nElements] : new vector<ZBosonInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEZBosonInfogR(void *p) {
      delete (static_cast<vector<ZBosonInfo>*>(p));
   }
   static void deleteArray_vectorlEZBosonInfogR(void *p) {
      delete [] (static_cast<vector<ZBosonInfo>*>(p));
   }
   static void destruct_vectorlEZBosonInfogR(void *p) {
      typedef vector<ZBosonInfo> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ZBosonInfo>

namespace ROOT {
   static TClass *vectorlESelectedMuongR_Dictionary();
   static void vectorlESelectedMuongR_TClassManip(TClass*);
   static void *new_vectorlESelectedMuongR(void *p = nullptr);
   static void *newArray_vectorlESelectedMuongR(Long_t size, void *p);
   static void delete_vectorlESelectedMuongR(void *p);
   static void deleteArray_vectorlESelectedMuongR(void *p);
   static void destruct_vectorlESelectedMuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<SelectedMuon>*)
   {
      vector<SelectedMuon> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<SelectedMuon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<SelectedMuon>", -2, "vector", 423,
                  typeid(vector<SelectedMuon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlESelectedMuongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<SelectedMuon>) );
      instance.SetNew(&new_vectorlESelectedMuongR);
      instance.SetNewArray(&newArray_vectorlESelectedMuongR);
      instance.SetDelete(&delete_vectorlESelectedMuongR);
      instance.SetDeleteArray(&deleteArray_vectorlESelectedMuongR);
      instance.SetDestructor(&destruct_vectorlESelectedMuongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<SelectedMuon> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<SelectedMuon>","std::vector<SelectedMuon, std::allocator<SelectedMuon> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<SelectedMuon>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlESelectedMuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<SelectedMuon>*>(nullptr))->GetClass();
      vectorlESelectedMuongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlESelectedMuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlESelectedMuongR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<SelectedMuon> : new vector<SelectedMuon>;
   }
   static void *newArray_vectorlESelectedMuongR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<SelectedMuon>[nElements] : new vector<SelectedMuon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlESelectedMuongR(void *p) {
      delete (static_cast<vector<SelectedMuon>*>(p));
   }
   static void deleteArray_vectorlESelectedMuongR(void *p) {
      delete [] (static_cast<vector<SelectedMuon>*>(p));
   }
   static void destruct_vectorlESelectedMuongR(void *p) {
      typedef vector<SelectedMuon> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<SelectedMuon>

namespace {
  void TriggerDictionaryInitialization_classes_Impl() {
    static const char* headers[] = {
"SelectedMuon.h",
"ZBosonInfo.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc7_amd64_gcc12/lcg/root/6.30.03-eeaa8ee64e8127bca194ba397d21d514/include/",
"/home/dndus0107/CMSSW_14_0_19_patch2/src/interface/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "classes dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$SelectedMuon.h")))  SelectedMuon;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate("$clingAutoload$ZBosonInfo.h")))  ZBosonInfo;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "classes dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "SelectedMuon.h"
#include "ZBosonInfo.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"SelectedMuon", payloadCode, "@",
"ZBosonInfo", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("classes",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_classes_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_classes_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_classes() {
  TriggerDictionaryInitialization_classes_Impl();
}
