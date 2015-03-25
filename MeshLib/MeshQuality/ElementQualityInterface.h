/**
 * \file   ElementQualityInterface.h
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the ElementQualityInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTQUALITYINTERFACE_H_
#define ELEMENTQUALITYINTERFACE_H_

#include <vector>

#include "BaseLib/Histogram.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshQuality/ElementQualityMetric.h"
#include "MeshLib/MeshQuality/EdgeRatioMetric.h"
#include "MeshLib/MeshQuality/ElementSizeMetric.h"
#include "MeshLib/MeshQuality/SizeDifferenceMetric.h"
#include "MeshLib/MeshQuality/AngleSkewMetric.h"
#include "MeshLib/MeshQuality/RadiusEdgeRatioMetric.h"

namespace MeshLib
{

/**
 * Interface class for handling mesh element quality metrics
 */
class ElementQualityInterface
{
public:
	/// Constructor
	ElementQualityInterface(MeshLib::Mesh const& mesh, MeshQualityType t)
	: _type(t), _mesh(mesh), _quality_tester(nullptr)
	{
		calculateElementQuality(_mesh, _type);
	}

	/// Destructor
	~ElementQualityInterface() 
	{ 
		delete _quality_tester; 
	}

	/// Returns the vector containing a quality measure for each element.
	std::vector<double> const& getQualityVector() const
	{ 
		return _quality_tester->getElementQuality(); 
	}

	/// Returns a histogram of the quality vector seperated into the given number of bins.
	/// If no number of bins is specified, one will be calculated based on the Sturges criterium.
	BaseLib::Histogram<double> getHistogram(std::size_t n_bins = 0) const
	{
		return _quality_tester->getHistogram(static_cast<size_t>(n_bins));
	}

	/// Writes a histogram of the quality vector to a specified file.
	int writeHistogram(std::string const& file_name, std::size_t n_bins = 0) const
	{
		if (_quality_tester == nullptr)
			return 1;

		BaseLib::Histogram<double> const histogram (_quality_tester->getHistogram(n_bins));
		histogram.write(file_name, _mesh.getName(), MeshQualityType2String(_type));
		return 0;
	}
	
private:
	/// Calculates the quality of each mesh element based on the specified metric
	void calculateElementQuality(MeshLib::Mesh const& mesh, MeshQualityType t)
	{
		if (t == MeshQualityType::EDGERATIO)
			_quality_tester = new MeshLib::EdgeRatioMetric(mesh);
		else if (t == MeshQualityType::ELEMENTSIZE)
			_quality_tester = new MeshLib::ElementSizeMetric(mesh);
		else if (t == MeshQualityType::SIZEDIFFERENCE)
			_quality_tester = new MeshLib::SizeDifferenceMetric(mesh);
		else if (t == MeshQualityType::EQUIANGLESKEW)
			_quality_tester = new MeshLib::AngleSkewMetric(mesh);
		else if (t == MeshQualityType::RADIUSEDGERATIO)
			_quality_tester = new MeshLib::RadiusEdgeRatioMetric(mesh);
		else
		{
			ERR("VtkVisPipeline::checkMeshQuality(): Unknown MeshQualityType.");
			return;
		}
		_quality_tester->calculateQuality();
	}

	MeshQualityType const _type;
	MeshLib::Mesh const& _mesh;
	MeshLib::ElementQualityMetric* _quality_tester;
};

}

#endif /* ELEMENTQUALITYINTERFACE_H_ */
