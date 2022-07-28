#include "Configurations.hpp"
#include "StringsUtilities.hpp"
#include "IOUtilities.hpp"
//#include "CsvExporter.hpp"
#include <set>

using namespace Eigen;

namespace Gedim
{
  // ***************************************************************************
  unordered_map<string, ConfigurationPropertySupportedTypes::SupportedTypes> Configurations::propertyTypes;
  unordered_map<string, void*> Configurations::properties;
  unsigned int Configurations::numberProperties = 0;
  typedef unsigned int uint;
  // ***************************************************************************
  void Configurations::RemoveProperty(const string& id,
                                      const ConfigurationPropertySupportedTypes::SupportedTypes& type)
  {
    switch (type)
    {
      case ConfigurationPropertySupportedTypes::Int:
        delete static_cast<const ConfigurationProperty<int>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::Char:
        delete static_cast<const ConfigurationProperty<char>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::String:
        delete static_cast<const ConfigurationProperty<string>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::Double:
        delete static_cast<const ConfigurationProperty<double>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::VectorInt:
        delete static_cast<const ConfigurationProperty<vector<int>>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::VectorDouble:
        delete static_cast<const ConfigurationProperty<vector<double>>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::Bool:
        delete static_cast<const ConfigurationProperty<bool>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::UInt:
        delete static_cast<const ConfigurationProperty<unsigned int>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::VectorUInt:
        delete static_cast<const ConfigurationProperty<vector<unsigned int>>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::SetUInt:
        delete static_cast<const ConfigurationProperty<set<unsigned int>>*>(properties[id]);
        break;
      case ConfigurationPropertySupportedTypes::Matrix:
        delete static_cast<const ConfigurationProperty<MatrixXd>*>(properties[id]);
        break;
      default:
        throw invalid_argument("Property '%s' cannot be deleted. Type not recognized");
        break;
    }

    properties.erase(id);
    numberProperties--;
  }
  // ***************************************************************************
  void Configurations::ConvertPropertyFromString(const string& id,
                                                 const string& type,
                                                 const string& value,
                                                 const string& description)
  {
    auto convertedType = ConfigurationPropertySupportedTypes::StringToType(type);

    switch (convertedType)
    {
      case ConfigurationPropertySupportedTypes::UInt:
      {
        AddProperty(id, uint(), description);

        unsigned int parsedValue = StringsUtilities::Parse<unsigned int>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::Int:
      {
        AddProperty(id, int(), description);

        int parsedValue = StringsUtilities::Parse<int>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::Char:
      {
        AddProperty(id, char(), description);

        char parsedValue = StringsUtilities::Parse<char>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::String:
      {
        AddProperty(id, string(), description);

        string parsedValue = StringsUtilities::Parse<string>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::Double:
      {
        AddProperty(id, double(), description);

        double parsedValue = StringsUtilities::Parse<double>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::VectorUInt:
      {
        AddProperty(id, vector<unsigned int>(), description);

        vector<unsigned int> parsedValue = StringsUtilities::Parse<vector<unsigned int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::VectorInt:
      {
        AddProperty(id, vector<int>(), description);

        vector<int> parsedValue = StringsUtilities::Parse<vector<int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::VectorDouble:
      {
        AddProperty(id, vector<double>(), description);

        vector<double> parsedValue = StringsUtilities::Parse<vector<double>>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::SetUInt:
      {
        AddProperty(id, set<unsigned int>(), description);

        set<unsigned int> parsedValue = StringsUtilities::Parse<set<unsigned int>>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      case ConfigurationPropertySupportedTypes::Bool:
      {
        AddProperty(id, bool(), description);

        bool parsedValue = StringsUtilities::Parse<bool>(value);
        SetPropertyValue(id, parsedValue);

        break;
      }
      default:
        throw invalid_argument("Property '" + id + "': type not supported");
    }
  }
  // ***************************************************************************
  void Configurations::AddProperty(const string& id,
                                   const char* value,
                                   const string& description)
  {
    RemoveProperty(id);

    properties.insert(make_pair(id, new ConfigurationProperty<string>()));
    propertyTypes.insert(make_pair(id,
                                   ConfigurationPropertySupportedTypes::String));
    numberProperties++;

    ConfigurationProperty<string>& property = GetProperty<string>(id);

    property.Id = id;
    property.Value = string(value);
    property.Description = description;
  }
  // ***************************************************************************
  void Configurations::SetPropertyValue(const string& id, const char* value)
  {
    if (!ExistsProperty(id))
      throw invalid_argument("");

    ConfigurationProperty<string>& property = GetProperty<string>(id);
    property.Value = string(value);
  }
  // ***************************************************************************
  void Configurations::InitializeFromCsv(const string& inputFile, const bool& hasHeader, const char& separator)
  {
    if (!Output::FileExists(inputFile))
      throw runtime_error("File '" + inputFile + "' not found");

    ifstream inFile;
    inFile.open(inputFile.c_str());

    if(!inFile.is_open())
      throw runtime_error("Cannot open file '" + inputFile + "'");

    int numberLines = 0;
    while (!inFile.eof())
    {
      string line;
      getline(inFile, line);

      // Skip Comment Line
      if(line[0] == '#')
        continue;

      if (numberLines == 0 && hasHeader)
      {
        numberLines++;
        continue;
      }

      vector<string> columns = StringsUtilities::Split(line, separator);

      if (columns.size() < 2)
        continue;

      switch (columns.size())
      {
        case 2:
        {
          const string& id = columns[0];
          const string& type = "string";
          const string& value = columns[1];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        case 3:
        {
          const string& id = columns[0];
          const string& type = columns[1];
          const string& value = columns[2];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        case 4:
        {
          const string& id = columns[0];
          const string& type = columns[1];
          const string& value = columns[2];
          const string& description = columns[3];

          ConvertPropertyFromString(id, type, value, description);
          break;
        }
        default:
          break;
      }
    }

    inFile.close();
  }
  // ***************************************************************************
  void Configurations::InitializeFromIni(const string& inputFile)
  {
    if (!Output::FileExists(inputFile))
      throw runtime_error("File '" + inputFile + "' not found");

    ifstream inFile;
    inFile.open(inputFile.c_str());

    if(!inFile.is_open())
      throw runtime_error("Cannot open file '" + inputFile + "'");

    // Find Property in File
    while (!inFile.eof())
    {
      string line;

      getline(inFile, line);

      // Skip Comment Line
      if(line[0] == '#')
        continue;

      vector<string> columns = StringsUtilities::Split(line, ' ');

      if (columns.size() < 2)
        continue;

      switch (columns.size())
      {
        case 2:
        {
          const string& id = columns[0];
          const string& type = "string";
          const string& value = columns[1];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        case 3:
        {
          const string& id = columns[0];
          const string& type = columns[1];
          const string& value = columns[2];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        default:
          break;
      }
    }

    inFile.close();


  }
  // ***************************************************************************
  void Configurations::Initialize(const int& argc,
                                  char **argv,
                                  const string& format)
  {
    char typeSeparator = StringsUtilities::FindSeparator(format, "id", "type");
    char valueSeparator = StringsUtilities::FindSeparator(format, "type", "value");

    for (int i = 0; i < argc; i++)
    {
      vector<string> strings = StringsUtilities::Split(argv[i], { typeSeparator, valueSeparator });

      switch (strings.size())
      {
        case 2:
        {
          const string& id = strings[0];
          const string& type = "string";
          const string& value = strings[1];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        case 3:
        {
          const string& id = strings[0];
          const string& type = strings[1];
          const string& value = strings[2];

          ConvertPropertyFromString(id, type, value);
          break;
        }
        default:
          break;
      }
    }


  }
  // ***************************************************************************
  void Configurations::ExportToCsv(const string& filePath,
                                   const bool& append,
                                   const char& separator)
  {
    throw runtime_error("Unimplemented method");

    // CsvExporter exporter(separator);

    // exporter.AddRowName("name");
    // exporter.AddRowName("type");
    // exporter.AddRowName("value");
    // exporter.AddRowName("description");

    // for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    // {
    //   const string& id = it->first;
    //   const ConfigurationPropertySupportedTypes::SupportedTypes type = it->second;

    //   switch (type)
    //   {
    //     case ConfigurationPropertySupportedTypes::Int:
    //       exporter.AddRow(GetProperty<int>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::Char:
    //       exporter.AddRow(GetProperty<char>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::String:
    //       exporter.AddRow(GetProperty<string>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::Double:
    //       exporter.AddRow(GetProperty<double>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::VectorInt:
    //       exporter.AddRow(GetProperty<vector<int>>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::VectorDouble:
    //       exporter.AddRow(GetProperty<vector<double>>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::Bool:
    //       exporter.AddRow(GetProperty<bool>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::UInt:
    //       exporter.AddRow(GetProperty<unsigned int>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::VectorUInt:
    //       exporter.AddRow(GetProperty<vector<unsigned int>>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::SetUInt:
    //       exporter.AddRow(GetProperty<set<unsigned int>>(id));
    //       break;
    //     case ConfigurationPropertySupportedTypes::Matrix:
    //       exporter.AddRow(GetProperty<MatrixXd>(id));
    //       break;
    //     default:
    //       Output::PrintWarningMessage("Property '%s' cannot be converted. Type not recognized", false);
    //       break;
    //   }
    // }

    // return exporter.Export(filePath, append);
  }
  // ***************************************************************************
  void Configurations::ExportToIni(const string& filePath,
                                   const bool& append,
                                   const string& section)
  {
    throw runtime_error("Unimplemented method");

    //    IniExporter exporter(section);

    //    for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    //    {
    //      const string& id = it->first;
    //      const ConfigurationPropertySupportedTypes::SupportedTypes type = it->second;

    //      switch (type)
    //      {
    //        case ConfigurationPropertySupportedTypes::Int:
    //          exporter.AddProperty(GetProperty<int>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::Char:
    //          exporter.AddProperty(GetProperty<char>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::String:
    //          exporter.AddProperty(GetProperty<string>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::Double:
    //          exporter.AddProperty(GetProperty<double>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::VectorInt:
    //          exporter.AddProperty(GetProperty<vector<int>>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::VectorDouble:
    //          exporter.AddProperty(GetProperty<vector<double>>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::Bool:
    //          exporter.AddProperty(GetProperty<bool>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::UInt:
    //          exporter.AddProperty(GetProperty<unsigned int>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::VectorUInt:
    //          exporter.AddProperty(GetProperty<vector<unsigned int>>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::SetUInt:
    //          exporter.AddProperty(GetProperty<set<unsigned int>>(id));
    //          break;
    //        case ConfigurationPropertySupportedTypes::Matrix:
    //          exporter.AddProperty(GetProperty<MatrixXd>(id));
    //          break;
    //        default:
    //          Output::PrintWarningMessage("Property '%s' cannot be converted. Type not recognized", false);
    //          break;
    //      }
    //    }

    //    return exporter.Export(filePath, append);
  }
  // ***************************************************************************
  void Configurations::Reset()
  {
    for (auto it = propertyTypes.begin(); it != propertyTypes.end(); it++)
    {
      const string& id = it->first;
      const ConfigurationPropertySupportedTypes::SupportedTypes& type = it->second;

      RemoveProperty(id, type);
    }

    propertyTypes.clear();
    numberProperties = 0;
  }
  // ***************************************************************************
  void Configurations::RemoveProperty(const string& id)
  {
    if (!ExistsProperty(id))
      return;

    RemoveProperty(id, GetPropertyType(id));
    propertyTypes.erase(id);
  }
  // ***************************************************************************
}
