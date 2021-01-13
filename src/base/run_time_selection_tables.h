#ifndef __RUN_TIME_SELECTION_TABLES_H__
#define __RUN_TIME_SELECTION_TABLES_H__

/*!
 * \brief 在所有 sub model 的类声明中声明 type name 以及函数
 * \param typeNameCString type name，类型应该是 const char*
 */
#define cfdemTypeName(typeNameCString)                       \
  static const char* cTypeName() { return typeNameCString; } \
  static const std::string typeName_;                        \
  virtual const std::string& typeName() const { return typeName_; }

/*!
 * \brief 在所有 base model 的类声明中声明 type name 以及函数
 * \param baseTypeNameCString base type name，类型应该是 const char*
 * \param defaultDerivedTypeNameCString default derived type name，类型应该是 const char*
 */
#define cfdemBaseTypeName(baseTypeNameCString, defaultDerivedTypeNameCString) \
  cfdemTypeName(baseTypeNameCString);                                         \
  static const char* cDefaultDerivedTypeName() { return defaultDerivedTypeNameCString; }

/*!
 * \brief 在所有 sub model 的类定义中初始化静态成员 typeName_
 *        无论是这个模型的基类，还是子类都需要声明
 * \param Type 指定当前类的类名
 */
#define cfdemDefineTypeName(Type) const std::string Type::typeName_(Type::cTypeName());

/*!
 * \brief 运行时选择器，这个宏只需要在每一个 sub model 的基类的声明中使用
 * \param AutoPtr 用于封装 sub model 指针的模板，通常就选择 Foam::autoPtr
 * \param BaseType 基类类型
 * \param argList 构造函数的形参
 * \param parList 构造函数的实参
 * @pre NewFunctionPtr: 定义一个函数指针类型，返回值类型为 autoPtr<baseType>，参数为 argList
 * @pre NewFunctionMap: 定义 unordered_map 的类型，key 类型为 std::string，value 类型为 NewFunctionPtr
 * @pre newFunctionMapPtr_: 声明 NewFunctionMap 类型静态指针
 * @pre constructNewFunctionMap: 创建 newFunctionMapPtr_ 的静态方法
 * @pre destroyNewFunctionMap: 析构 newFunctionMapPtr_ 的静态方法
 * @pre NewFunctionAdder: 在子类的实现中创建 NewFunctionAdder<子类> 的一个实例
 *                        则这个实例在构造方法中，会将 NewFunctionAdder<子类>::New 这个静态方法插入到 map 中，
 *                        而这个静态方法就是实际用来创建子类对象的方法
 */
#define cfdemDeclareRunTimeSelection(AutoPtr, BaseType, argList, parList)                                           \
  typedef AutoPtr<BaseType> NewFunctionReturnType;                                                                  \
                                                                                                                    \
  typedef NewFunctionReturnType(*NewFunctionPtr) argList;                                                           \
                                                                                                                    \
  typedef std::unordered_map<std::string, NewFunctionPtr> NewFunctionMap;                                           \
                                                                                                                    \
  static NewFunctionMap* newFunctionMapPtr_;                                                                        \
                                                                                                                    \
  static void constructNewFunctionMap();                                                                            \
                                                                                                                    \
  static void destroyNewFunctionMap();                                                                              \
                                                                                                                    \
  template <typename BaseType##Type>                                                                                \
  class NewFunctionAdder {                                                                                          \
   public:                                                                                                          \
    static AutoPtr<BaseType> NewFunction argList {                                                                  \
      if (false == std::is_base_of<BaseType, BaseType##Type>::value) {                                              \
        std::cerr << "New error: " << BaseType##Type::typeName_ << " is not convertible to " << #BaseType           \
                  << std::endl;                                                                                     \
        error::safePrintStack(std::cerr);                                                                           \
      }                                                                                                             \
      return AutoPtr<BaseType>(new BaseType##Type parList);                                                         \
    }                                                                                                               \
    NewFunctionAdder(const std::string& typeName = BaseType##Type::typeName_) {                                     \
      constructNewFunctionMap();                                                                                    \
      std::cout << "Add New function from " << typeName << " to " << #BaseType << "'s map" << std::endl;            \
      if (newFunctionMapPtr_->find(typeName) == newFunctionMapPtr_->end()) {                                        \
        newFunctionMapPtr_->emplace(typeName, NewFunction);                                                         \
      } else {                                                                                                      \
        std::cerr << "Duplicate entry: " << typeName << " has been added to runtime selection map in " << #BaseType \
                  << std::endl;                                                                                     \
        error::safePrintStack(std::cerr);                                                                           \
      }                                                                                                             \
    }                                                                                                               \
    ~NewFunctionAdder() { destroyNewFunctionMap(); }                                                                \
  };

/*!
 * \brief 定义 BaseType 中的静态成员变量 newFunctionMapPtr_
 * \param BaseType 指定基类
 */
#define cfdemDefineNewFunctionMap(BaseType) BaseType::NewFunctionMap* BaseType::newFunctionMapPtr_ = nullptr;

/*!
 * \brief 定义 BaseType 中的静态方法 constructNewFunctionMap()，用于创建对象 map 对象
 * \param BaseType 指定基类
 */
#define cfdemDefineConstructNewFunctionMap(BaseType)               \
  void BaseType::constructNewFunctionMap() {                       \
    static bool constructed = false;                               \
    if (!constructed) {                                            \
      constructed = true;                                          \
      BaseType::newFunctionMapPtr_ = new BaseType::NewFunctionMap; \
    }                                                              \
  }

/*!
 * \brief 定义 BaseType 中的静态方法 newRunTimeSelectionConstructorMap()，用于创建对象 map 对象
 * \param BaseType 指定基类
 */
#define cfdemDefineDestroyNewFunctionMap(BaseType) \
  void BaseType::destroyNewFunctionMap() {         \
    if (BaseType::newFunctionMapPtr_ != nullptr) { \
      delete BaseType::newFunctionMapPtr_;         \
      BaseType::newFunctionMapPtr_ = nullptr;      \
    }                                              \
  }

/*!
 * \brief 用于在基类中定义 New 这个静态方法，而这个方法返回的就是运行期间动态选择的某一个子类对象
 * \param AutoPtr 用于封装 sub model 指针的模板，通常就选择 Foam::autoPtr
 * \param BaseType 指定基类
 * \param NewArgList New 这个静态方法的形参
 * \param paramTypeName NewArgList 中用于指定模型名称的参数
 * \param paramConstructorList 构造函数的实参
 */
#define cfdmeDefineBaseTypeNewWithTypeName(AutoPtr, BaseType, NewArgList, paramTypeName, paramConstructorList) \
  AutoPtr<BaseType> BaseType::New NewArgList {                                                                 \
    std::string name(paramTypeName);                                                                           \
    const char* pName = name.c_str();                                                                          \
    Info << "\nSelecting " << #BaseType << ": " << pName << endl;                                              \
    if (BaseType::newFunctionMapPtr_->find(pName) == BaseType::newFunctionMapPtr_->end()) {                    \
      FatalError << #BaseType << "::New" << #NewArgList << ": "                                                \
                 << "unknow " << #BaseType << " type " << pName << ", constructor not in map" << endl          \
                 << "Valid " << #BaseType << " types are: [";                                                  \
      for (const std::pair<std::string, BaseType::NewFunctionPtr>& item : *(BaseType::newFunctionMapPtr_)) {   \
        FatalError << item.first << ", ";                                                                      \
      }                                                                                                        \
      FatalError << "]" << endl << abort(FatalError);                                                          \
    }                                                                                                          \
    BaseType::NewFunctionPtr selectNewFunction = (BaseType::newFunctionMapPtr_)->operator[](pName);            \
    Info << "Construct " << pName << " model..." << endl;                                                      \
    return selectNewFunction paramConstructorList;                                                             \
  }

/*!
 * \brief 用于在基类中定义 New 这个静态方法，而这个方法返回的就是运行期间动态选择的某一个子类对象
 * \param AutoPtr 用于封装 sub model 指针的模板，通常就选择 Foam::autoPtr
 * \param BaseType 指定基类
 * \param NewArgList New 这个静态方法的形参
 * \param paramDict NewArgList 中 const dicttionary& 类型的实参
 * \param paramConstructorList 构造函数的实参
 */
#define cfdmeDefineBaseTypeNew(AutoPtr, BaseType, NewArgList, paramDict, paramConstructorList)                        \
  AutoPtr<BaseType> BaseType::New NewArgList {                                                                        \
    std::string name =                                                                                                \
        paramDict.found(#BaseType) ? word(paramDict.lookup(#BaseType)).c_str() : BaseType::cDefaultDerivedTypeName(); \
    const char* pName = name.c_str();                                                                                 \
    Info << "\nSelecting " << #BaseType << ": " << pName << endl;                                                     \
    if (BaseType::newFunctionMapPtr_->find(pName) == BaseType::newFunctionMapPtr_->end()) {                           \
      FatalError << #BaseType << "::New" << #NewArgList << ": "                                                       \
                 << "unknow " << #BaseType << " type " << pName << ", constructor not in map" << endl                 \
                 << "Valid " << #BaseType << " types are: [";                                                         \
      for (const std::pair<std::string, BaseType::NewFunctionPtr>& item : *(BaseType::newFunctionMapPtr_)) {          \
        Info << item.first << ", ";                                                                                   \
      }                                                                                                               \
      Info << "]" << endl << abort(FatalError);                                                                       \
    }                                                                                                                 \
    BaseType::NewFunctionPtr selectNewFunction = (BaseType::newFunctionMapPtr_)->operator[](pName);                   \
    Info << "Construct " << pName << " model..." << endl;                                                             \
    return selectNewFunction paramConstructorList;                                                                    \
  }

#define cfdemDefineNewFunctionAdder(BaseType, DerivedType) \
  static BaseType::NewFunctionAdder<DerivedType> DerivedType##NewFunctionAdder_;

#define cfdemCreateNewFunctionAdder(BaseType, DerivedType) \
  BaseType::NewFunctionAdder<DerivedType> DerivedType::DerivedType##NewFunctionAdder_;

#endif  // __RUN_TIME_SELECTION_TABLES_H__
