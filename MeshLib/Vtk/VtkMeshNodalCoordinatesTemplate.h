/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  VtkMeshNodalCoordinatesTemplate is a adapter for node coordinates of
 *         OGS meshes to VTK unstructured grids.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vtkMappedDataArray.h>
#include <vtkObjectFactory.h>  // for vtkStandardNewMacro
#include <vtkVersion.h>

#if VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1
#include <vtkTypeTemplate.h>   // For templated vtkObject API
#endif

namespace MeshLib
{
    class Node;
}

namespace MeshLib
{
template <class Scalar>
class VtkMeshNodalCoordinatesTemplate :
#if !(VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1)
    public vtkMappedDataArray<Scalar>
#else
    public vtkTypeTemplate<VtkMeshNodalCoordinatesTemplate<Scalar>,
                           vtkMappedDataArray<Scalar>>
#endif // vtk version
{
public:
#if !(VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1)
    vtkTemplateTypeMacro(VtkMeshNodalCoordinatesTemplate<Scalar>,
                         vtkMappedDataArray<Scalar>);
#else
    vtkMappedDataArrayNewInstanceMacro(VtkMeshNodalCoordinatesTemplate<Scalar>);
#endif // vtk version
    static VtkMeshNodalCoordinatesTemplate *New();
    void PrintSelf(std::ostream& os, vtkIndent indent) override;

    /// Pass the nodes from OGS mesh
    void SetNodes(std::vector<MeshLib::Node*> const & nodes);

    // Reimplemented virtuals -- see superclasses for descriptions
    void Initialize() override;
    void GetTuples(vtkIdList *ptIds, vtkAbstractArray *output) override;
    void GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output) override;
    void Squeeze() override;
    vtkArrayIterator *NewIterator() override;
    vtkIdType LookupValue(vtkVariant value) override;
    void LookupValue(vtkVariant value, vtkIdList *ids) override;
    vtkVariant GetVariantValue(vtkIdType idx) override;
    void ClearLookup() override;
    double* GetTuple(vtkIdType i) override;
    void GetTuple(vtkIdType i, double *tuple) override;
    vtkIdType LookupTypedValue(Scalar value) override;
    void LookupTypedValue(Scalar value, vtkIdList *ids) override;
    Scalar& GetValueReference(vtkIdType idx) override;

    // This container is read only -- this method does nothing but print a
    // warning.
    int Allocate(vtkIdType sz, vtkIdType ext) override;
    int Resize(vtkIdType numTuples) override;
    void SetNumberOfTuples(vtkIdType number) override;
    void SetTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
    void SetTuple(vtkIdType i, const float *source) override;
    void SetTuple(vtkIdType i, const double *source) override;
    void InsertTuple(vtkIdType i, vtkIdType j, vtkAbstractArray *source) override;
    void InsertTuple(vtkIdType i, const float *source) override;
    void InsertTuple(vtkIdType i, const double *source) override;
    void InsertTuples(vtkIdList *dstIds, vtkIdList *srcIds,
                      vtkAbstractArray *source) override;
    void InsertTuples(vtkIdType /*unused*/, vtkIdType /*unused*/,
                      vtkIdType /*unused*/,
                      vtkAbstractArray* /*unused*/) override;
    vtkIdType InsertNextTuple(vtkIdType j, vtkAbstractArray *source) override;
    vtkIdType InsertNextTuple(const float *source) override;
    vtkIdType InsertNextTuple(const double *source) override;
    void InsertVariantValue(vtkIdType idx, vtkVariant value) override;
    void DeepCopy(vtkAbstractArray *aa) override;
    void DeepCopy(vtkDataArray *da) override;
    void InterpolateTuple(vtkIdType i, vtkIdList *ptIndices,
                          vtkAbstractArray* source,  double* weights) override;
    void InterpolateTuple(vtkIdType i, vtkIdType id1, vtkAbstractArray *source1,
                          vtkIdType id2, vtkAbstractArray *source2, double t) override;
    void SetVariantValue(vtkIdType idx, vtkVariant value) override;
    void RemoveTuple(vtkIdType id) override;
    void RemoveFirstTuple() override;
    void RemoveLastTuple() override;
    void SetValue(vtkIdType idx, Scalar value) override;
    vtkIdType InsertNextValue(Scalar v) override;
    void InsertValue(vtkIdType idx, Scalar v) override;

#if !(VTK_MAJOR_VERSION < 7 || VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION < 1)
    Scalar& GetValueReference(vtkIdType idx) const;
    Scalar GetValue(vtkIdType idx) const override;
    void GetTypedTuple(vtkIdType tupleId, Scalar* t) const override;
    void SetTypedTuple(vtkIdType i, const Scalar* t) override;
    void InsertTypedTuple(vtkIdType i, const Scalar* t) override;
    vtkIdType InsertNextTypedTuple(const Scalar* t) override;
#else
    Scalar GetValue(vtkIdType idx) override;
    void GetTupleValue(vtkIdType idx, Scalar* t) override;
    void SetTupleValue(vtkIdType i, const Scalar* t) override;
    void InsertTupleValue(vtkIdType i, const Scalar* t) override;
    vtkIdType InsertNextTupleValue(const Scalar* t) override;
#endif  // vtk version

    VtkMeshNodalCoordinatesTemplate(const VtkMeshNodalCoordinatesTemplate&) =
        delete;
    void operator=(const VtkMeshNodalCoordinatesTemplate&) = delete;

protected:
    VtkMeshNodalCoordinatesTemplate();
    ~VtkMeshNodalCoordinatesTemplate() override;

    const std::vector<MeshLib::Node*>* _nodes{nullptr};

private:
    vtkIdType Lookup(const Scalar &val, vtkIdType startIndex);
    double* TempDoubleArray{nullptr};
};

}  // namespace MeshLib

#include "VtkMeshNodalCoordinatesTemplate-impl.h"
