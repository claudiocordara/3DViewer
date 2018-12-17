#ifndef SCENE3D_H
#define SCENE3D_H
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/Subdivision_method_3.h>
#include <CGAL/assertions.h>
#include <CGAL/exceptions.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h> 
#include <list>

#include <vector>
#include <QtGui> 
#include <QtWidgets>
#include <QtOpenGL/QGLWidget>


class MainWindow;

// define the most used CGAL types
typedef CGAL::Simple_cartesian<double>									Kernel;
typedef Kernel::Point_3													Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>	Polyhedron;
typedef Polyhedron::Facet_iterator										Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator					Halfedge_facet_circulator;
typedef Polyhedron::Halfedge_around_facet_const_circulator				Halfedge_facet_const_circulator;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor				vertex_descriptor;


#define			SEGGRAPH_NOTMARKED				0
#define			SEGGRAPH_VISITED				1

class SegmentGraph {
public:

	class SegmentNode {
	public:
		SegmentNode();
		SegmentNode(Polyhedron *_seg);
		~SegmentNode();

		Polyhedron* segment() const;
		void segment(Polyhedron* _seg);
		int childCount() const;
		SegmentNode* getChildAt(int _index) const;
		void appendChild(SegmentNode* _node);
		void insertChildAfter(SegmentNode* _node, int _index);
		SegmentNode* removeChildAt(int _index);
		void clearChildList();
		int getFlag() const;
		void setFlag(int _flag);
	private:
		int mFlag;
		Polyhedron *pSegment;
		std::vector<SegmentNode*> mChildList;
	};

	SegmentGraph();
	~SegmentGraph();

	void addSegment(Polyhedron *_seg);
	SegmentNode* getNodeWithSegment(Polyhedron *_seg) const;
	void addConnection(Polyhedron *_seg1, Polyhedron *_seg2);
	void addChild(Polyhedron *_father, Polyhedron *_child);
	void addChild(SegmentNode *_father, SegmentNode *_child);
	void addHRoot(Polyhedron *_seg);
	void addHRoot(SegmentNode *_node);
	void clear();
	int nodeCount() const;
	SegmentNode* getNodeAt(int _index) const;
	std::vector<SegmentNode*> getConnectedNodes(SegmentNode *_node) const;
	int hRootCount() const;
	SegmentNode* getHRootAt(int _index) const;
private:
	std::vector<SegmentNode*> mNodeList;
	std::vector<std::pair<SegmentNode*,SegmentNode*> > mEdgeList;
	std::vector<SegmentNode*> mHRootList;
};



// Scene3D class to 3D objects visualization using Qt
class Scene3D : public QGLWidget
{
private:
	GLfloat xRot;				// the rotation angle of X axis
	GLfloat yRot;				// the rotation angle of Y axis
	GLfloat zRot;				// the rotation angle of Z axis
	GLfloat zTra;				// translation by Z axis
	GLfloat nSca;				// scale of image

	int vmesh;					// the number of mesh vertices
	int fmesh;					// the number of mesh facets
	int vpart;					// the number of parts vertices
	int fpart;					// the number of parts facets
	int vskel;					// the number of skeleton vertices
	int eskel;					// the number of skeleton edges
	int maxgroup;				// the number of skeleton's unjoined groups

	GLfloat asp;				// aspect ratio to set the maximum distance from mesh surface to origin equal to one
	// shifts to set the center of mesh on the point of origin
	GLfloat cx;					// shift by X axis
	GLfloat cy;					// shifh by Y axis
	GLfloat cz;					// shift by Z axis
	QPoint ptrMousePosition;	// the last saved mouse position

	GLfloat *MeshVertex;		// array of solid's vertices

	GLfloat *MeshColor;			// array of solid's colors of vertices (for facets)
	GLfloat *MeshColorSeg;		// array of solid's colors of vertices (for facets)
	GLfloat *MeshColorSDF;		// array of solid's colors of vertices (for facets)

	GLuint *MeshIndex;			// array of solid's facets
	GLfloat *EdgeColor;			// array of solid's colors of vertices (for edges)
	GLuint *EdgeIndex;			// array of solid's edges
	GLfloat *SkelVertex;		// array of skeleton's vertices
	GLfloat *SkelSegmColors;	// array of skeleton's random colors of segments
	GLfloat *SkelColor;			// array of skeleton's colors of vertices (for edges)
	GLuint *SkelIndex;			// array of skeleton's edges
	GLint *SkelSegVert;			// array of skeleton's relations vertex<->segment
	GLuint *SkelSegEdge;		// array of skeleton's edges (doubled)
	GLuint *SkelMapMesh;		// array of relations skeleton<->solid vertices
	GLfloat *PartVertex;		// array of part's vertices
	GLfloat *PartColor;			// array of part's colors of vertices (for facets)
	GLuint *PartIndex;			// array of part's facets
	GLfloat *PEdgeColor;		// array of part's colors of vertices (for edges)
	GLuint *PEdgeIndex;			// array of part's edges
	std::vector<int> SegGroups;	// array of skeleton's unjoined parts
	SegmentGraph SegGraph;

	MainWindow *mParentWnd;

	void scale_plus();
	void scale_minus();
	void rotate_up();
	void rotate_down();
	void rotate_left();
	void rotate_right();
	void translate_down();
	void translate_up();
	void defaultScene();

	void drawAxis();
	void getAsp();
	void getMesh();
	void addMesh(Polyhedron P);
	void drawSkeleton();
	void drawWireframe();
	void drawFacets();

protected:
	void initializeGL();
	void resizeGL(int nWidth, int nHeight);
	void paintGL();
	void mousePressEvent(QMouseEvent* pe);
	void mouseMoveEvent(QMouseEvent* pe);
	void mouseReleaseEvent(QMouseEvent* pe);
	void wheelEvent(QWheelEvent* pe);

public:
	int showElem;
	std::vector<Polyhedron> tmesh;

	std::vector<Polyhedron> segMesh;
	
	Scene3D(MainWindow * parent);
	void load(Polyhedron mesh);
	void add(Polyhedron mesh);
	void updateMesh();
	void keyPressEvent(QKeyEvent* pe)
	{
		// set the actions of keyboard keys
		switch (pe->key())
		{
		case Qt::Key_Plus: scale_plus(); break;
		case Qt::Key_Equal: scale_plus(); break;
		case Qt::Key_Minus: scale_minus(); break;
		case Qt::Key_Up: rotate_up(); break;
		case Qt::Key_Down: rotate_down(); break;
		case Qt::Key_Left: rotate_left(); break;
		case Qt::Key_Right: rotate_right(); break;
		case Qt::Key_Z: translate_down(); break;
		case Qt::Key_X: translate_up(); break;
		case Qt::Key_Space: defaultScene(); break;
		case Qt::Key_Escape: this->close(); break;
		}
		updateGL();
	}
	int getSkeleton();
	void switchColors();
	void getParts();

	void newPart();
	bool splitPartsBySkeleton();
	bool splitPartBySDF();
	bool splitPartBySkeletonAndSDF();


	
	void groupsToOff();
	void drawParts();
	void drawPWireframe();
	int testSegmentationBySDF();
	int testSegmentationBySkeletonAndSDF();
	int testPolyedraDecomposition();

	int computeSegmentHierarchy();
	int updateSegmentationTree();
};


typedef boost::graph_traits<Polyhedron>::face_descriptor             face_descriptor;

// Property map associating a facet with an integer as id to an
// element in a vector stored internally
template<class ValueType>
struct Facet_with_id_pmap : public boost::put_get_helper<ValueType&, Facet_with_id_pmap<ValueType> > {
	typedef face_descriptor key_type;
	typedef ValueType value_type;
	typedef value_type& reference;
	typedef boost::lvalue_property_map_tag category;
	Facet_with_id_pmap(
		std::vector<ValueType>& internal_vector
	) : internal_vector(internal_vector) { }
	reference operator[](key_type key) const
	{
		return internal_vector[key->id()];
	}
private:
	std::vector<ValueType>& internal_vector;
};


#endif
