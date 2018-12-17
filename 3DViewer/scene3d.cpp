#include <fstream>
#include <vector>
#include <QtGui> 
#include <QtWidgets>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include "scene3D.h"
#include "functions.h"
#include "mainwindow.h"

// define the most used CGAL types
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>			Skeletonization;
typedef Skeletonization::Skeleton										Skeleton;
typedef Skeleton::vertex_descriptor										Skeleton_vertex;
typedef Skeleton::edge_descriptor										Skeleton_edge;


typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef K::Segment_3 Segment;

// Initiation of Scene3D object
Scene3D::Scene3D(MainWindow* parent) : QGLWidget(parent)
{
	mParentWnd = parent;
	// default values of parameters
	xRot = -90; yRot = 0; zRot = 0; zTra = 0; nSca = 1;	vmesh = 0; vskel = 0; vpart = 0;
}

// Initiation of OpenGL
void Scene3D::initializeGL()
{
	// set the background color
	qglClearColor(Qt::white);
	// to hide the covered deeper objects
	glEnable(GL_DEPTH_TEST);
	// disable the shade
	glShadeModel(GL_FLAT);
	// to use the blended colors
	glEnable(GL_BLEND);
	// enable the transparency
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// calculate the aspect ratio
	getAsp();
	// calculate the mesh arrays
	getMesh();
	// calculate the parts arrays
	getParts();
	// to use the arrays of vertices for drawing
	glEnableClientState(GL_VERTEX_ARRAY);
	// to use the arrays of colors for drawing
	glEnableClientState(GL_COLOR_ARRAY);
}

// Recalculate the scene parameters after the window will resized
void Scene3D::resizeGL(int nWidth, int nHeight)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// choose the most problem direction
	GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;
	if (nWidth >= nHeight)
		// change the matrix using horizontal ratio
		glOrtho(-1.0 / ratio, 1.0 / ratio, -1.0, 1.0, -10.0, 1.0);
	else
		// change the matrix using vertical ratio
		glOrtho(-1.0, 1.0, -1.0*ratio, 1.0*ratio, -10.0, 1.0);
	// transform to window coordinates
	glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

// Draw the scene
void Scene3D::paintGL()
{
	// set the initial parameters
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	// reset the transformation matrix to identity
	glLoadIdentity();
	// apply the scale
	glScalef(nSca, nSca, nSca);
	// apply the translation
	glTranslatef(0.0f, zTra, 0.0f);
	// apply the rotations
	glRotatef(xRot, 1.0f, 0.0f, 0.0f);
	glRotatef(yRot, 0.0f, 1.0f, 0.0f);
	glRotatef(zRot, 0.0f, 0.0f, 1.0f);

	// draw the elements using the 'elements visibility' variable
	if (showElem & shAxis) drawAxis();
	if (showElem & shSkeleton) drawSkeleton();
	if (showElem & shWireframe) drawWireframe();
	if ((showElem & shWireframe) && (showElem & shParts)) drawPWireframe();
	if (showElem & shFacets) drawFacets();
	if (showElem & shParts) drawParts();
}

// Set the initial position of actions doing by mouse
void Scene3D::mousePressEvent(QMouseEvent* pe)
{
	// save the mouse position
	ptrMousePosition = pe->pos();
}

// Event don't used; it was created for consistent purpose
void Scene3D::mouseReleaseEvent(QMouseEvent* pe)
{
	//
}

// Process rotation by mouse
void Scene3D::mouseMoveEvent(QMouseEvent* pe)
{
	// calculate rotation by X axis
	xRot += 180 / nSca * (GLfloat)(pe->y() - ptrMousePosition.y()) / height();
	// calculate rotation by Z axis
	zRot += 180 / nSca * (GLfloat)(pe->x() - ptrMousePosition.x()) / width();
	// save the mouse position
	ptrMousePosition = pe->pos();
	// draw the scene
	updateGL();
}

// Process zoom by mouse
void Scene3D::wheelEvent(QWheelEvent* pe)
{
	// choose the scale method
	if ((pe->delta()) > 0) scale_plus(); else if ((pe->delta()) < 0) scale_minus();
	// draw the scene
	updateGL();
}

// The controls actions of the scene
// zoom
void Scene3D::scale_plus() { nSca = nSca * 1.1; }
void Scene3D::scale_minus() { nSca = nSca/1.1; }
// rotation
void Scene3D::rotate_up() { xRot += 1.0f; }
void Scene3D::rotate_down() { xRot -= 1.0f; }
void Scene3D::rotate_left() { zRot += 1.0f; }
void Scene3D::rotate_right() { zRot -= 1.0f; }
// translation
void Scene3D::translate_down() { zTra -= 0.05f; }
void Scene3D::translate_up() { zTra += 0.05f; }
// reset to default
void Scene3D::defaultScene() { xRot=-90; yRot=0; zRot=0; zTra=0; nSca=1; }

// Draw the axis
void Scene3D::drawAxis() 
{
	// set the line width
	glLineWidth(3.0f);
	// set the mode to 'lines'
	glBegin(GL_LINES);
	
	// the X axis
	glColor4f(1.00f, 0.00f, 0.00f, 1.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);
	glVertex3f(-1.0f, 0.0f, 0.0f);
	// the Y axis
	glColor4f(0.00f, 1.00f, 0.00f, 1.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, -1.0f, 0.0f);
	// the Z axis
	glColor4f(0.00f, 0.00f, 1.00f, 1.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, -1.0f);
	
	// end of the 'lines' mode
	glEnd();
}

// Calculate the aspect ratio of given mesh
void Scene3D::getAsp()
{
	// initiate the ratio and center point
	asp = 1;
	cx = 0;
	cy = 0;
	cz = 0;

	// initiate the total number of vertices of all parts
	vmesh = 0;
	// if we have no parts return
	int lmesh = (int)tmesh.size();
	if (lmesh == 0) return;

	// loop over the parts
	for (int p = 0; p < tmesh.size(); p++)
		// add the number of vertices
		vmesh += (int)tmesh[p].size_of_vertices();
	// if we have no vertices return
	if (vmesh == 0) return;

	// initiate the maximum and minimum values of X, Y and Z
	int i = 0;
	float mx, my, mz;
	float Mx, My, Mz;
	// loop over the parts
	for (int p = 0; p < lmesh; p++) {
		// loop over the vertices
		for (auto it = tmesh[p].vertices_begin(); it != tmesh[p].vertices_end(); it++) {
			// check and replace the maximum and minimum values of X, Y and Z
			auto p = it->point();
			if (i == 0) { mx = p.x(); Mx = p.x(); }
			if (mx > p.x()) mx = p.x();
			if (Mx < p.x()) Mx = p.x();
			if (i == 0) { my = p.y(); My = p.y(); }
			if (my > p.y()) my = p.y();
			if (My < p.y()) My = p.y();
			if (i == 0) { mz = p.z(); Mz = p.z(); }
			if (mz > p.z()) mz = p.z();
			if (Mz < p.z()) Mz = p.z();
			i++;
		}
	}

	// define the middle point
	cx = (mx + Mx) / 2;
	cy = (my + My) / 2;
	cz = (mz + Mz) / 2;

	// calculate the lowest ratio
	if (Mx != mx) asp = 1 / (Mx - mx);
	if (My != my) { if (1 / (My - my) < asp) asp = 1 / (My - my); }
	if (Mz != mz) { if (1 / (Mz - mz) < asp) asp = 1 / (Mz - mz); }
}

// Calculate the arrays of mesh to use OpenGL (vertices, colors, facets, edges)
void Scene3D::getMesh()
{
	// check we have a correct mesh
	vmesh = 0;
	int lmesh = (int)tmesh.size() - 1;
	if (lmesh < 0) return;
	if (tmesh[lmesh].size_of_vertices() == 0) return;

	// get the number of facets and vertices
	vmesh = (int)tmesh[lmesh].size_of_vertices();
	fmesh = (int)tmesh[lmesh].size_of_facets();

	// resize the arrays
	if (MeshVertex != NULL) {
		delete MeshVertex;
	}
	MeshVertex = new GLfloat[3 * vmesh];
	if (MeshColor != NULL) {
		delete MeshColor;
	}
	MeshColor = new GLfloat[4 * vmesh];
	if (MeshIndex != NULL) {
		delete MeshIndex;
	}
	MeshIndex = new GLuint[3 * fmesh];
	if (EdgeColor != NULL) {
		delete EdgeColor;
	}
	EdgeColor = new GLfloat[4 * vmesh];
	if (EdgeIndex != NULL) {
		delete EdgeIndex;
	}
	EdgeIndex = new GLuint[6 * fmesh];
	if (SkelMapMesh != NULL) {
		delete SkelMapMesh;
	}
	SkelMapMesh = new GLuint[vmesh];
	if (MeshColorSDF != NULL) {
		delete MeshColorSDF;
	}
	MeshColorSDF = new GLfloat[4 * vmesh];
	if (MeshColorSeg != NULL) {
		delete MeshColorSeg;
	}
	MeshColorSeg = new GLfloat[4 * vmesh];

	// loop over vertices
	int i = 0;
	for (auto it = tmesh[lmesh].vertices_begin(); it != tmesh[lmesh].vertices_end(); it++, i++) {
		auto p = it->point();
		it->id() = i;
		// add vertex coordinates to array
		MeshVertex[3 * i] = (p.x() - cx) * asp;
		MeshVertex[3 * i + 1] = (p.y() - cy) * asp;
		MeshVertex[3 * i + 2] = (p.z() - cz) * asp;

		// add vertex colors to array of facets colors
		MeshColor[4 * i] = 0.5f;
		MeshColor[4 * i + 1] = 0.7f;
		MeshColor[4 * i + 2] = 0.5f;
		MeshColor[4 * i + 3] = 0.3f;

		MeshColorSDF[4 * i] = 0.5f;
		MeshColorSDF[4 * i + 1] = 0.7f;
		MeshColorSDF[4 * i + 2] = 0.5f;
		MeshColorSDF[4 * i + 3] = 0.3f;

		// add vertex color to array of edges colors
		EdgeColor[4 * i] = 0.25f;
		EdgeColor[4 * i + 1] = 0.0f;
		EdgeColor[4 * i + 2] = 0.0f;
		EdgeColor[4 * i + 3] = 1.0f;
	}

	// loop over facets
	i = 0;
	for (Facet_iterator fi = tmesh[lmesh].facets_begin(); fi != tmesh[lmesh].facets_end(); ++fi, i++) {
		// get the facet's vertices enumerator
		Halfedge_facet_circulator fj = fi->facet_begin();
		CGAL_assertion(CGAL::circulator_size(fj) >= 3);
		// loop over facet's vertices
		int j = 0;
		do
			// add vertices index to array
			MeshIndex[3 * i + j++] = std::distance(tmesh[lmesh].vertices_begin(), fj->vertex());
		while (++fj != fi->facet_begin());
	}
	// loop over facets
	for (i = 0; i < fmesh; i++) {
		// add the edges to array
		EdgeIndex[6 * i] = MeshIndex[3 * i];
		EdgeIndex[6 * i + 1] = MeshIndex[3 * i + 1];
		EdgeIndex[6 * i + 2] = MeshIndex[3 * i + 1];
		EdgeIndex[6 * i + 3] = MeshIndex[3 * i + 2];
		EdgeIndex[6 * i + 4] = MeshIndex[3 * i + 2];
		EdgeIndex[6 * i + 5] = MeshIndex[3 * i];
	}
}

// Calculate the arrays of mesh to use OpenGL (vertices, colors, facets, edges)
void Scene3D::getParts()
{
	SegmentGraph::SegmentNode* selectedSeg = mParentWnd->getSelectedSegment();
	if (selectedSeg != NULL && selectedSeg->segment() != NULL) {
		// check do we have the correct parts
		//vpart = 0;
		//int lmesh = (int)tmesh.size() - 1;
		//if (lmesh < 1) return;
		//for (int p = 0; p < lmesh; p++)
		//	vpart += (int)tmesh[p].size_of_vertices();
		//if (vpart == 0) return;

		//// create the common mesh from the parts
		//Polyhedron mesh = mergePoly(tmesh, 0, lmesh);

		//// get the number of facets and vertices
		//vpart = (int)mesh.size_of_vertices();
		//fpart = (int)mesh.size_of_facets();
		Polyhedron *poly = selectedSeg->segment();
		vpart = (int)poly->size_of_vertices();
		fpart = (int)poly->size_of_facets();

		// resize the arrays
		if (PartVertex != NULL) {
			delete PartVertex;
		}
		PartVertex = new GLfloat[3 * vpart];
		if (PartColor != NULL) {
			delete PartColor;
		}
		PartColor = new GLfloat[4 * vpart];
		if (PartIndex != NULL) {
			delete PartIndex;
		}
		PartIndex = new GLuint[3 * fpart];
		if (PEdgeColor != NULL) {
			delete PEdgeColor;
		}
		PEdgeColor = new GLfloat[4 * vpart];
		if (PEdgeIndex != NULL) {
			delete PEdgeIndex;
		}
		PEdgeIndex = new GLuint[6 * fpart];

		// loop over vertices
		int i = 0;
		for (auto it = poly->vertices_begin(); it != poly->vertices_end(); it++, i++) {
			auto p = it->point();
			it->id() = i;

			// add vertex coordinates to array
			PartVertex[3 * i] = (p.x() - cx) * asp;
			PartVertex[3 * i + 1] = (p.y() - cy) * asp;
			PartVertex[3 * i + 2] = (p.z() - cz) * asp;

			// add vertex colors to array of facet's colors
			PartColor[4 * i] = 1.0f;
			PartColor[4 * i + 1] = 0.0f;
			PartColor[4 * i + 2] = 0.0f;
			PartColor[4 * i + 3] = 1.0f;

			// add vertex colors to array of edge's colors
			PEdgeColor[4 * i] = 1.0f;
			PEdgeColor[4 * i + 1] = 0.0f;
			PEdgeColor[4 * i + 2] = 0.0f;
			PEdgeColor[4 * i + 3] = 1.0f;
		}

		// loop over facets
		i = 0;
		for (Facet_iterator fi = poly->facets_begin(); fi != poly->facets_end(); ++fi, i++) {
			// get the facet's vertices enumerator
			Halfedge_facet_circulator fj = fi->facet_begin();
			CGAL_assertion(CGAL::circulator_size(fj) >= 3);
			// loop over facet's vertices
			int j = 0;
			do
				// add vertex indices to array
				PartIndex[3 * i + j++] = std::distance(poly->vertices_begin(), fj->vertex());
			while (++fj != fi->facet_begin());
		}
		// loop over facets
		for (i = 0; i < fpart; i++) {
			// add the edges to array
			PEdgeIndex[6 * i] = PartIndex[3 * i];
			PEdgeIndex[6 * i + 1] = PartIndex[3 * i + 1];
			PEdgeIndex[6 * i + 2] = PartIndex[3 * i + 1];
			PEdgeIndex[6 * i + 3] = PartIndex[3 * i + 2];
			PEdgeIndex[6 * i + 4] = PartIndex[3 * i + 2];
			PEdgeIndex[6 * i + 5] = PartIndex[3 * i];
		}
	}
	else {
		// check do we have the correct parts
		vpart = 0;
		int lmesh = (int)tmesh.size() - 1;
		if (lmesh < 1) return;
		for (int p = 0; p < lmesh; p++)
			vpart += (int)tmesh[p].size_of_vertices();
		if (vpart == 0) return;

		// create the common mesh from the parts
		Polyhedron mesh = mergePoly(tmesh, 0, lmesh);

		// get the number of facets and vertices
		vpart = (int)mesh.size_of_vertices();
		fpart = (int)mesh.size_of_facets();

		// resize the arrays
		PartVertex = new GLfloat[3 * vpart];
		PartColor = new GLfloat[4 * vpart];
		PartIndex = new GLuint[3 * fpart];
		PEdgeColor = new GLfloat[4 * vpart];
		PEdgeIndex = new GLuint[6 * fpart];

		// loop over vertices
		int i = 0;
		for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); it++, i++) {
			auto p = it->point();
			it->id() = i;

			// add vertex coordinates to array
			PartVertex[3 * i] = (p.x() - cx) * asp;
			PartVertex[3 * i + 1] = (p.y() - cy) * asp;
			PartVertex[3 * i + 2] = (p.z() - cz) * asp;

			// add vertex colors to array of facet's colors
			PartColor[4 * i] = 0.7f;
			PartColor[4 * i + 1] = 0.5f;
			PartColor[4 * i + 2] = 0.5f;
			PartColor[4 * i + 3] = 0.3f;

			// add vertex colors to array of edge's colors
			PEdgeColor[4 * i] = 0.0f;
			PEdgeColor[4 * i + 1] = 0.0f;
			PEdgeColor[4 * i + 2] = 0.25f;
			PEdgeColor[4 * i + 3] = 1.0f;
		}

		// loop over facets
		i = 0;
		for (Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); ++fi, i++) {
			// get the facet's vertices enumerator
			Halfedge_facet_circulator fj = fi->facet_begin();
			CGAL_assertion(CGAL::circulator_size(fj) >= 3);
			// loop over facet's vertices
			int j = 0;
			do
				// add vertex indices to array
				PartIndex[3 * i + j++] = std::distance(mesh.vertices_begin(), fj->vertex());
			while (++fj != fi->facet_begin());
		}
		// loop over facets
		for (i = 0; i < fpart; i++) {
			// add the edges to array
			PEdgeIndex[6 * i] = PartIndex[3 * i];
			PEdgeIndex[6 * i + 1] = PartIndex[3 * i + 1];
			PEdgeIndex[6 * i + 2] = PartIndex[3 * i + 1];
			PEdgeIndex[6 * i + 3] = PartIndex[3 * i + 2];
			PEdgeIndex[6 * i + 4] = PartIndex[3 * i + 2];
			PEdgeIndex[6 * i + 5] = PartIndex[3 * i];
		}
	}

}

// Create the CGAL skeleton of the main segment and calculate its arrays
int Scene3D::getSkeleton()
{
	// check do we have the correct main part
	int lmesh = (int)tmesh.size() - 1;
	if (lmesh < 0) return 0;
	if (vmesh == 0) return 0;
	if (tmesh[lmesh].size_of_vertices() == 0) return 0;
	
	// calculate skeleton using CGAL
	Skeleton skeleton;
	CGAL::extract_mean_curvature_flow_skeleton(tmesh[lmesh], skeleton);

	// get the number of facets and vertices
	vskel = (int)boost::num_vertices(skeleton);
	eskel = (int)boost::num_edges(skeleton);

	// resize the arrays
	SkelVertex = new GLfloat[3 * vskel];
	SkelColor = new GLfloat[4 * vskel];
	SkelIndex = new GLuint[2 * eskel];
	SkelSegVert = new GLint[vskel];
	SkelSegEdge = new GLuint[4 * eskel];

	// loop over skeleton's vertices
	int i = 0;
	BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton)) {
		// loop over joined solid's vertices
		BOOST_FOREACH(vertex_descriptor vd, skeleton[v].vertices)
			// add relation to array
			SkelMapMesh[vd->id()] = (GLuint)v;
		
		// add vertex coordinates to array
		SkelVertex[3 * i] = (skeleton[v].point.x() - cx) * asp;
		SkelVertex[3 * i + 1] = (skeleton[v].point.y() - cy) * asp;
		SkelVertex[3 * i + 2] = (skeleton[v].point.z() - cz) * asp;
		
		SkelSegVert[i] = -1;
		i++;
	}

	// loop over skeleton's edges
	i = 0;
	BOOST_FOREACH(Skeleton_edge e, edges(skeleton)) {
		// get the indices of the start and end points
		GLuint s = (GLuint)e.m_source;
		GLuint t = (GLuint)e.m_target;

		// add vertex indices to arrays
		SkelIndex[2 * i] = s;
		SkelIndex[2 * i + 1] = t;
		SkelSegEdge[4 * i] = s;
		SkelSegEdge[4 * i + 1] = t;
		SkelSegEdge[4 * i + 2] = t;
		SkelSegEdge[4 * i + 3] = s;
		i++;
	}

	// prepare the segmentation of skeleton
	SegGroups = SkelSeg(vskel, eskel, SkelSegEdge, SkelSegVert);
	// get the number of unjoined skeleton's parts
	int maxSeg = (int)SegGroups.size();
	maxgroup = SegGroups[maxSeg - 1];

	// resize the segments colors array
	SkelSegmColors = new GLfloat[3 * maxSeg];
	// loop over segments
	for (i = 0; i < maxSeg; i++) {
		int mainCol = i % 3;
		// fill the colors with random (except of blue)
		//SkelSegmColors[3 * i] = 0.05f + 0.01f * (qrand() % 15);
		//SkelSegmColors[3 * i + 1] = 0.05f + 0.01f * (qrand() % 15);
		//SkelSegmColors[3 * i + 2] = 0.05f + 0.01f * (qrand() % 15);
		//SkelSegmColors[3 * i + mainCol] += 0.75f;
		SkelSegmColors[3 * i] = (double)(qrand() % 255) / 255.0;
		SkelSegmColors[3 * i + 1] = (double)(qrand() % 255) / 255.0;
		SkelSegmColors[3 * i + 2] = (double)(qrand() % 255) / 255.0;
	}

	switchColors();
	/*
	// loop over skeleton's vertices
	for (i = 0; i < vskel; i++) {
		
		// blue color by default for skeleton
		GLfloat R = 0.0f, G = 0.0f, B = 0.75f;
		if (showElem & shColors) {
			R = SkelSegmColors[3 * SkelSegVert[i]];
			G = SkelSegmColors[3 * SkelSegVert[i] + 1];
			B = SkelSegmColors[3 * SkelSegVert[i] + 2];
		}
		// add vertex colors to array of edge's colors
		SkelColor[4 * i] = R;
		SkelColor[4 * i + 1] = G;
		SkelColor[4 * i + 2] = B;
		SkelColor[4 * i + 3] = 1.00f;
	}*/

	// loop over unjoined parts; try to find the parts which can be joined
	while (maxgroup > 0)
	{
		// initiate the shortest edge's index and length
		int bestEdge = -1;
		double bestDist = 100;

		// loop over solid's edges
		for (i = 0; i < 3 * fmesh; i++) {
			// get the unjoined part index of start and end of the edge
			int grS = SegGroups[SkelSegVert[SkelMapMesh[EdgeIndex[2 * i]]]];
			int grT = SegGroups[SkelSegVert[SkelMapMesh[EdgeIndex[2 * i + 1]]]];
			
			// if they are different
			if (grS != grT) {
				// calculate the length of the edge
				double dx = MeshVertex[3 * EdgeIndex[2 * i]] - MeshVertex[3 * EdgeIndex[2 * i + 1]];
				double dy = MeshVertex[3 * EdgeIndex[2 * i] + 1] - MeshVertex[3 * EdgeIndex[2 * i + 1] + 1];
				double dz = MeshVertex[3 * EdgeIndex[2 * i] + 2] - MeshVertex[3 * EdgeIndex[2 * i + 1] + 2];
				double curDist = dx * dx + dy * dy + dz * dz;
				
				// if this edge is shorter than the previous best
				if (curDist < bestDist) {
					// change the index and length of the shortest edge
					bestDist = curDist;
					bestEdge = i;
				}
			}
		}

		// we have no possibility to join some pair of parts
		if (bestEdge == -1) break;

		// get the unjoined part index of start and end of the best edge
		int grL = SegGroups[SkelSegVert[SkelMapMesh[EdgeIndex[2 * bestEdge]]]];
		int grH = SegGroups[SkelSegVert[SkelMapMesh[EdgeIndex[2 * bestEdge + 1]]]];
		// use the higher index
		if (grL > grH) grH = grL;

		// add the new sceleton's edge
		eskel++;
		SkelSegEdge = (GLuint *)realloc(SkelSegEdge, 4 * eskel * sizeof(GLuint));
		SkelSegEdge[4 * eskel - 4] = SkelMapMesh[EdgeIndex[2 * bestEdge]];
		SkelSegEdge[4 * eskel - 3] = SkelMapMesh[EdgeIndex[2 * bestEdge + 1]];
		SkelSegEdge[4 * eskel - 2] = SkelMapMesh[EdgeIndex[2 * bestEdge + 1]];
		SkelSegEdge[4 * eskel - 1] = SkelMapMesh[EdgeIndex[2 * bestEdge]];
		
		// loop over unjoined parts
		for (i = 0; i < maxSeg; i++)
			// shift all the indices higher than cutted
			if (SegGroups[i] >= grH) SegGroups[i]--;
		maxgroup--;
	}
	
	// we still have unjoined parts
	if (maxgroup > 0)
	{
		// create the error message
		std::stringstream ss;
		ss << "The model consist of " << maxgroup + 1 << " solids.";
		// show the message
		QMessageBox msgBox;
		msgBox.setText(QString::fromStdString(ss.str()));
		msgBox.setInformativeText("Do you want to save them as separate OFF files?");
		msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Cancel);
		msgBox.setDefaultButton(QMessageBox::Save);
		int ret = msgBox.exec();
		
		// if user choose the 'Save' then save unjoined parts as separate OFF files
		if (ret == QMessageBox::Save) groupsToOff();
	}
	
	// return the number of unjoined parts: zero is normal return
	return maxgroup;
}

// Switch segments colors between random and fixed
void Scene3D::switchColors()
{
	// check do the skeleton exist
	//if (vskel == 0) return;

	// loop over skeleton's vertices
	for (int i = 0; i < vskel; i++) {

		// blue color by default for skeleton
		GLfloat R = 0.0f, G = 0.0f, B = 0.75f;
		// set the different colors to segments if checked
		if (showElem & shColors) {
			R = SkelSegmColors[3 * SkelSegVert[i]];
			G = SkelSegmColors[3 * SkelSegVert[i] + 1];
			B = SkelSegmColors[3 * SkelSegVert[i] + 2];
		}
		// add vertex colors to array of edge's colors
		SkelColor[4 * i] = R;
		SkelColor[4 * i + 1] = G;
		SkelColor[4 * i + 2] = B;
		SkelColor[4 * i + 3] = 1.00f;
	}

	// loop over mesh vertices
	for (int i = 0; i < vmesh; i++) {
		// light-green transparent color by default
		GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
		if ((showElem & shColors) && vskel > 0) {
			R = SkelSegmColors[3 * SkelSegVert[SkelMapMesh[i]]];
			G = SkelSegmColors[3 * SkelSegVert[SkelMapMesh[i]] + 1];
			B = SkelSegmColors[3 * SkelSegVert[SkelMapMesh[i]] + 2];
			T = 1.0f;
		}
		else if (showElem & segColors) {
			R = MeshColorSeg[4 * i];
			G = MeshColorSeg[4 * i + 1];
			B = MeshColorSeg[4 * i + 2];
			T = MeshColorSeg[4 * i + 3];
		}
		else if (showElem & sdfColors) {
			R = MeshColorSDF[4 * i];
			G = MeshColorSDF[4 * i + 1];
			B = MeshColorSDF[4 * i + 2];
			T = MeshColorSDF[4 * i + 3];
		}

		// add vertex colors to array of facets colors
		MeshColor[4 * i] = R;
		MeshColor[4 * i + 1] = G;
		MeshColor[4 * i + 2] = B;
		MeshColor[4 * i + 3] = T;
	}
}



bool Scene3D::splitPartsBySkeleton() {
	std::vector<std::vector<double> > partsCoords;
	std::vector<std::vector<int> > pTris;			// part's facets
	std::vector<std::vector<int> > pVert;			// part's vertices (relation to old indices)
	bool ret = false;

	if (SkelSegEdge == NULL || SkelSegVert == NULL || SkelVertex == NULL) {
		getSkeleton();
	}

	if (SkelSegEdge != NULL && SkelSegVert != NULL && SkelVertex != NULL) {
		// segs[i][0] = Id of the i-th segment the segment will contains all skeleton vertices with SkelSegVert[vertex index] == Id
		// segs[i][1] = Index of the end point of the segment, it is a skeleton vertex index
		// segs[i][2] = Number of vertices in the segment
		// segs[i][3] = Index of the start point of the segment, it is a skeleton vertex index, is -1 if it is an outer segment	
		std::vector<std::array<int32_t, 4>> segs = GetAllSeg(vskel, eskel, SkelSegEdge, SkelSegVert);

		//std::vector<int> mVert;			// main vertices (relation to old indices)
		if (segs.empty() == false && segs[0][0] != -1) {
			ret = true;

			partsCoords.resize(segs.size());
			pVert.resize(segs.size());
			pTris.resize(segs.size());
			for (int i = 0; i < (int)segs.size(); i++) {
				for (int j = 0; j < vmesh + 1; j++) {
					// initiate the arrays of indices relations 
					pVert[i].push_back(-1);
					//mVert.push_back(-1);
				}
			}

			for (int i = 0; i < (int)segs.size(); i++) {
				if (segs[i][1] != -1) {
					for (int v = 0; v < 3; v++) {
						partsCoords[i].push_back(SkelVertex[3 * segs[i][1] + v]);
					}
				}
				if (segs[i][3] != -1) {
					for (int v = 0; v < 3; v++) {
						partsCoords[i].push_back(SkelVertex[3 * segs[i][3] + v]);
					}
				}
			}

			std::vector<int> pvert;
			for (int i = 0 ; i < (int)segs.size() ; i++) {
				pvert.push_back(1);
			}
			//int mvert = 1;


			for (int j = 0; j < fmesh && ret == true; j++) {
				int segIdVertex0 = SkelSegVert[SkelMapMesh[MeshIndex[3 * j]]];
				int segIdVertex1 = SkelSegVert[SkelMapMesh[MeshIndex[3 * j + 1]]];
				int segIdVertex2 = SkelSegVert[SkelMapMesh[MeshIndex[3 * j + 2]]];
				int segIdFace = -1;
				if (segIdVertex0 == segIdVertex1 || segIdVertex0 == segIdVertex2) {
					segIdFace = segIdVertex0;
				}
				else if (segIdVertex1 == segIdVertex0 || segIdVertex1 == segIdVertex2) {
					segIdFace = segIdVertex1;
				}
				else if (segIdVertex2 == segIdVertex0 || segIdVertex2 == segIdVertex1) {
					segIdFace = segIdVertex2;
				}

				if (segIdFace != -1) {
					int segNdx = -1;
					for (int j = 0; j < (int)segs.size() && segNdx == -1; j++) {
						if (segIdFace == segs[j][0]) {
							segNdx = j;
						}
					}
					if (segNdx != -1) {
						// loop over facet's vertices
						for (int v = 0; v < 3; v++) {
							// get the index of current vertex
							int vn = MeshIndex[3 * j + v];
							// it was not saved in array yet
							if (pVert[segNdx][vn] == -1)
							{
								// add the relation between solid's and part's indices
								pVert[segNdx][vn] = pvert[segNdx]++;// pvert++;
								// add the coordinates to part's array


								partsCoords[segNdx].push_back(MeshVertex[3 * vn]);
								partsCoords[segNdx].push_back(MeshVertex[3 * vn + 1]);
								partsCoords[segNdx].push_back(MeshVertex[3 * vn + 2]);
							}
							// add the new index to part's array
							pTris[segNdx].push_back(pVert[segNdx][vn]);
						}
					}
					else {
						ret = false;
					}
				}
				else {
					ret = false;
				}
			}


			try {
				segMesh.clear();
				for (int i = 0; i < (int)segs.size(); i++) {
					Polyhedron part = arrToMesh(partsCoords[i], pTris[i]);
					segMesh.push_back(part);
				}
				computeSegmentHierarchy();
				updateSegmentationTree();
				//for (int i = 0; i < (int)segMesh.size(); i++) {
				//	// create the serial part of filename
				//	std::stringstream sc;
				//	sc << "F:\\Sviluppo\\3DViewer\\Data\\part" << i << ".off";
				//	meshToOff((char*)(sc.str().c_str()), segMesh[i]);
				//}
			}
			catch (const CGAL::Assertion_exception e) { // the extracted vertices and facets don't form a proper mesh
				// create the error message
				QMessageBox msg;
				msg.setText("The part extraction was unsuccessful, please reopen the model and try again");
				msg.setStandardButtons(QMessageBox::Ok);
				msg.setDefaultButton(QMessageBox::Ok);
				msg.exec();
			}
		}
	}

	return ret;
}


bool Scene3D::splitPartBySDF() {
	bool ret = false;
	if (tmesh.size() == 1) {
		Polyhedron mesh = tmesh[0];
		// Segmentation from sdf values
		// create a property-map for SDF values
		typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
		Facet_double_map internal_sdf_map;
		boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);
		// compute SDF values using default parameters for number of rays, and cone angle

		CGAL::sdf_values(mesh,
			sdf_property_map,
			2.0 / 3.0 * CGAL_PI,	// double cone_angle = 2.0 / 3.0 * CGAL_PI
			25,						// std::size_t number_of_rays = 25
			true					// bool postprocess = true
		);


		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			double value = sdf_property_map[facet_it];
			if (0 <= value && value <= 1.0 / 8.0) {
				R = 0.0;
				G = 0.0;
				B = 4.0 * value + .5; // .5 - 1 // b = 1/2
			}
			else if (1.0 / 8.0 < value && value <= 3.0 / 8.0) {
				R = 0.0;
				G = 4.0 * value - .5; // 0 - 1 // b = - 1/2
				B = 1.0; // small fix
			}
			else if (3.0 / 8.0 < value && value <= 5.0 / 8.0) {
				R = 4.0 * value - 1.5; // 0 - 1 // b = - 3/2
				G = 1.0;
				B = -4.0 * value + 2.5; // 1 - 0 // b = 5/2
			}
			else if (5.0 / 8.0 < value && value <= 7.0 / 8.0) {
				R = 1.0;
				G = -4.0 * value + 3.5; // 1 - 0 // b = 7/2
				B = 0.0;
			}
			else if (7.0 / 8.0 < value && value <= 1.0) {
				R = -4.0 * value + 4.5; // 1 - .5 // b = 9/2
				G = 0.0;
				B = 0.0;
			}
			else {    // should never happen - value > 1
				R = .5;
				G = 0;
				B = 0;
			}
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				MeshColorSDF[4 * vIndex] = R;
				MeshColorSDF[4 * vIndex + 1] = G;
				MeshColorSDF[4 * vIndex + 2] = B;
				MeshColorSDF[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());
		}



		// create a property-map for segment-ids
		typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
		Facet_int_map internal_segment_map;
		boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

		// Segment the mesh using default parameters for number of levels, and smoothing lambda
		// Any other scalar values can be used instead of using SDF values computed using the CGAL function
		std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh,
			sdf_property_map,
			segment_property_map,
			5,							// std::size_t number_of_clusters = 5
			0.26,						// double smoothing_lambda = 0.26
			false						// bool output_cluster_ids = false
		);

		std::vector<GLfloat> segmColorsArray;// [3 * number_of_segments];
		segmColorsArray.resize(3 * number_of_segments);
		// loop over segments
		for (int i = 0; i < number_of_segments; i++) {
			int mainCol = i % 3;
			// fill the colors with random (except of blue)
			segmColorsArray[3 * i] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 1] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 2] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + mainCol] += 0.75f;
		}


		// print segment-ids

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin() ; facet_it != mesh.facets_end(); ++facet_it) {
			//std::cout << segment_property_map[facet_it] << " ";

			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			// set the different colors to segments if checked
			R = segmColorsArray[3 * segment_property_map[facet_it]];
			G = segmColorsArray[3 * segment_property_map[facet_it] + 1];
			B = segmColorsArray[3 * segment_property_map[facet_it] + 2];
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			//std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSeg[4 * vIndex] = R;
				MeshColorSeg[4 * vIndex + 1] = G;
				MeshColorSeg[4 * vIndex + 2] = B;
				MeshColorSeg[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());

		}
		//std::cout << std::endl;

		if (number_of_segments > 0) {
			std::vector<std::vector<double> > partsCoords;
			std::vector<std::vector<int> > pTris;			// part's facets
			std::vector<std::vector<int> > pVert;			// part's vertices (relation to old indices)

			partsCoords.resize(number_of_segments);
			pVert.resize(number_of_segments);
			pTris.resize(number_of_segments);
			for (int i = 0; i < number_of_segments; i++) {
				for (int j = 0; j < vmesh + 1; j++) {
					// initiate the arrays of indices relations 
					pVert[i].push_back(-1);
				}
			}

			std::vector<int> pvert;
			for (int i = 0; i < number_of_segments; i++) {
				pvert.push_back(0);
			}

			for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
				int faceSegmentIndex = (int)segment_property_map[facet_it];

				Halfedge_facet_circulator j = facet_it->facet_begin();
				do {
					
					int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
					if (pVert[faceSegmentIndex][vIndex] == -1) {
						pVert[faceSegmentIndex][vIndex] = pvert[faceSegmentIndex]++;// pvert++;

						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).x());
						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).y());
						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).z());

					}
					pTris[faceSegmentIndex].push_back(pVert[faceSegmentIndex][vIndex]);
				}
				while (++j != facet_it->facet_begin());
			}

			try {

				//std::vector<double> testCoord;
				//std::vector<int> testTris;
				//testCoord.push_back(0); testCoord.push_back(0); testCoord.push_back(0);
				//testCoord.push_back(10); testCoord.push_back(0); testCoord.push_back(0);
				//testCoord.push_back(5); testCoord.push_back(10); testCoord.push_back(0);
				//testTris.push_back(0);
				//testTris.push_back(1);
				//testTris.push_back(2);
				//Polyhedron testPart = arrToMesh(testCoord, testTris);
				//meshToOff("F:\\Sviluppo\\3DViewer\\Data\\testPart.off", testPart);

				segMesh.clear();
				for (int i = 0; i < number_of_segments; i++) {
					Polyhedron part = arrToMesh(partsCoords[i], pTris[i]);
					segMesh.push_back(part);
				}

				computeSegmentHierarchy();
				updateSegmentationTree();

				//for (int i = 0; i < (int)segMesh.size(); i++) {
				//	// create the serial part of filename
				//	std::stringstream sc;
				//	sc << "F:\\Sviluppo\\3DViewer\\Data\\part" << i << ".off";
				//	meshToOff((char*)(sc.str().c_str()), segMesh[i]);
				//}
			}
			catch (const CGAL::Assertion_exception e) { // the extracted vertices and facets don't form a proper mesh
				// create the error message
				QMessageBox msg;
				msg.setText("The part extraction was unsuccessful, please reopen the model and try again");
				msg.setStandardButtons(QMessageBox::Ok);
				msg.setDefaultButton(QMessageBox::Ok);
				msg.exec();
			}
		}
		ret = true;
	}

	return ret;
}


bool Scene3D::splitPartBySkeletonAndSDF() {
	bool ret = false;
	if (tmesh.size() == 1) {
		Polyhedron mesh = tmesh[0];

		// extract the skeleton
		Skeleton skeleton;
		CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
		// init the polyhedron simplex indices
		CGAL::set_halfedgeds_items_id(mesh);
		//for each input vertex compute its distance to the skeleton
		std::vector<double> distances(num_vertices(mesh));
		BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
		{
			const Point& skel_pt = skeleton[v].point;
			BOOST_FOREACH(vertex_descriptor mesh_v, skeleton[v].vertices)
			{
				const Point& mesh_pt = mesh_v->point();
				distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
			}
		}

		// create a property-map for sdf values
		std::vector<double> sdf_values(num_faces(mesh));
		Facet_with_id_pmap<double> sdf_property_map(sdf_values);
		// compute sdf values with skeleton
		BOOST_FOREACH(face_descriptor f, faces(mesh))
		{
			double dist = 0;
			BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, mesh), mesh))
				dist += distances[target(hd, mesh)->id()];
			sdf_property_map[f] = dist / 3.;
		}
		// post-process the sdf values
		CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
			//std::cout << sdf_property_map[facet_it] << " ";
			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			double value = sdf_property_map[facet_it];
			if (0 <= value && value <= 1.0 / 8.0) {
				R = 0.0;
				G = 0.0;
				B = 4.0 * value + .5; // .5 - 1 // b = 1/2
			}
			else if (1.0 / 8.0 < value && value <= 3.0 / 8.0) {
				R = 0.0;
				G = 4.0 * value - .5; // 0 - 1 // b = - 1/2
				B = 1.0; // small fix
			}
			else if (3.0 / 8.0 < value && value <= 5.0 / 8.0) {
				R = 4.0 * value - 1.5; // 0 - 1 // b = - 3/2
				G = 1.0;
				B = -4.0 * value + 2.5; // 1 - 0 // b = 5/2
			}
			else if (5.0 / 8.0 < value && value <= 7.0 / 8.0) {
				R = 1.0;
				G = -4.0 * value + 3.5; // 1 - 0 // b = 7/2
				B = 0.0;
			}
			else if (7.0 / 8.0 < value && value <= 1.0) {
				R = -4.0 * value + 4.5; // 1 - .5 // b = 9/2
				G = 0.0;
				B = 0.0;
			}
			else {    // should never happen - value > 1
				R = .5;
				G = 0;
				B = 0;
			}
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSDF[4 * vIndex] = R;
				MeshColorSDF[4 * vIndex + 1] = G;
				MeshColorSDF[4 * vIndex + 2] = B;
				MeshColorSDF[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());
		}


		// create a property-map for segment-ids (it is an adaptor for this case)
		std::vector<std::size_t> segment_ids(num_faces(mesh));
		Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
		// Segment the mesh using default parameters for number of levels, and smoothing lambda
		// Any other scalar values can be used instead of using SDF values computed using the CGAL function
		std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh,
			sdf_property_map,
			segment_property_map,
			5,							// std::size_t number_of_clusters = 5
			0.26,						// double smoothing_lambda = 0.26
			false						// bool output_cluster_ids = false
		);



		std::vector<GLfloat> segmColorsArray;// [3 * number_of_segments];
		segmColorsArray.resize(3 * number_of_segments);
		// loop over segments
		for (int i = 0; i < number_of_segments; i++) {
			int mainCol = i % 3;
			// fill the colors with random (except of blue)
			segmColorsArray[3 * i] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 1] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 2] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + mainCol] += 0.75f;
		}


		// print segment-ids

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
			//std::cout << segment_property_map[facet_it] << " ";

			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			// set the different colors to segments if checked
			R = segmColorsArray[3 * segment_property_map[facet_it]];
			G = segmColorsArray[3 * segment_property_map[facet_it] + 1];
			B = segmColorsArray[3 * segment_property_map[facet_it] + 2];
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			//std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSeg[4 * vIndex] = R;
				MeshColorSeg[4 * vIndex + 1] = G;
				MeshColorSeg[4 * vIndex + 2] = B;
				MeshColorSeg[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());

		}
		//std::cout << std::endl;

		if (number_of_segments > 0) {
			std::vector<std::vector<double> > partsCoords;
			std::vector<std::vector<int> > pTris;			// part's facets
			std::vector<std::vector<int> > pVert;			// part's vertices (relation to old indices)

			partsCoords.resize(number_of_segments);
			pVert.resize(number_of_segments);
			pTris.resize(number_of_segments);
			for (int i = 0; i < number_of_segments; i++) {
				for (int j = 0; j < vmesh + 1; j++) {
					// initiate the arrays of indices relations 
					pVert[i].push_back(-1);
				}
			}

			std::vector<int> pvert;
			for (int i = 0; i < number_of_segments; i++) {
				pvert.push_back(0);
			}

			for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
				int faceSegmentIndex = (int)segment_property_map[facet_it];

				Halfedge_facet_circulator j = facet_it->facet_begin();
				do {

					int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
					if (pVert[faceSegmentIndex][vIndex] == -1) {
						pVert[faceSegmentIndex][vIndex] = pvert[faceSegmentIndex]++;// pvert++;

						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).x());
						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).y());
						partsCoords[faceSegmentIndex].push_back((j->vertex()->point()).z());

					}
					pTris[faceSegmentIndex].push_back(pVert[faceSegmentIndex][vIndex]);
				} while (++j != facet_it->facet_begin());
			}

			try {
				segMesh.clear();
				for (int i = 0; i < number_of_segments; i++) {
					Polyhedron part = arrToMesh(partsCoords[i], pTris[i]);
					segMesh.push_back(part);
				}

				computeSegmentHierarchy();
				updateSegmentationTree();

				//for (int i = 0; i < (int)segMesh.size(); i++) {
				//	// create the serial part of filename
				//	std::stringstream sc;
				//	sc << "F:\\Sviluppo\\3DViewer\\Data\\part" << i << ".off";
				//	meshToOff((char*)(sc.str().c_str()), segMesh[i]);
				//}
			}
			catch (const CGAL::Assertion_exception e) { // the extracted vertices and facets don't form a proper mesh
				// create the error message
				QMessageBox msg;
				msg.setText("The part extraction was unsuccessful, please reopen the model and try again");
				msg.setStandardButtons(QMessageBox::Ok);
				msg.setDefaultButton(QMessageBox::Ok);
				msg.exec();
			}
		}
		ret = true;
	}

	return ret;
}


// Extract the new segment from solid
void Scene3D::newPart()
{
	// get the longest 'output branch' of skeleton
	std::array<int32_t, 3> segs = GetSeg(vskel, eskel, SkelSegEdge, SkelSegVert);
	// check do we have correct branch
	if (segs[0] >= 0 && segs[2] > 2) {
		
		// initiate the arrays of new extracted part and of the rest of solid
		std::vector<int> pVert;			// part's vertices (relation to old indices)
		std::vector<double> pCoords;	// part's vertices coordinates
		std::vector<int> pTris;			// part's facets
		std::vector<int> mVert;			// main vertices (relation to old indices)
		std::vector<double> mCoords;	// main vertices coordinates
		std::vector<int> mTris;			// main facets

		// loop over solid's vertices  + new vertex from skeleton
		for (int j = 0; j < vmesh + 1; j++) {
			// initiate the arrays of indices relations 
			pVert.push_back(-1);
			mVert.push_back(-1);
		}
		// create the new vertex from last vertex of skeleton's branch which will be used to cut solid
		for (int v = 0; v < 3; v++) {
			// add it to both part's and main arrays
			pCoords.push_back(SkelVertex[3 * segs[1] + v]);
			mCoords.push_back(SkelVertex[3 * segs[1] + v]);
		}
		// initiate the counters of number of vertices (part's and main)
		int pvert = 1;
		int mvert = 1;

		// loop over solid's facets
		for (int j = 0; j < fmesh; j++)
		{
			// initiate the checkers
			int isSeg = 0;		// number of vertices which is of new part
			int notSeg = -1;	// index of vertex which is not of the new part

			// loop over facet's vertices
			for (int v = 0; v < 3; v++) {
				// set the checkers
				if (SkelSegVert[SkelMapMesh[MeshIndex[3 * j + v]]] == segs[0]) isSeg++;
				else notSeg = v;
			}
			
			// all the vertices is of the new part
			if (isSeg == 3) {
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in array yet
					if (pVert[vn] == -1)
					{
						// add the relation between solid's and part's indices
						pVert[vn] = pvert++;
						// add the coordinates to part's array
						pCoords.push_back(MeshVertex[3 * vn]);
						pCoords.push_back(MeshVertex[3 * vn + 1]);
						pCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to part's array
					pTris.push_back(pVert[vn]);
				}
			}
			
			// the facet is on the border between part and main
			else if (isSeg == 2)
			{
				// save the existing faset to main
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in array yet
					if (mVert[vn] == -1) {
						// add the relation between solid's and main indices
						mVert[vn] = mvert++;
						// add the coordinates to main array
						mCoords.push_back(MeshVertex[3 * vn]);
						mCoords.push_back(MeshVertex[3 * vn + 1]);
						mCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to main array
					mTris.push_back(mVert[vn]);
				}

				// create new facets to main and part
				// initiate the temporary facet
				std::array<int, 3> vr;
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];

					// change the vertex which is not of new part to skeleton's vertex
					if (notSeg == v) {
						// into part's array
						pTris.push_back(0);
						// into temporary facet (for main array)
						vr[v] = 0;
					}
					// the vertices which is on new part will stay
					else {
						// it was not saved in part's array yet
						if (pVert[vn] == -1) {
							// add the relation between solid's and part's indices
							pVert[vn] = pvert++;
							// add the coordinates to part's array
							pCoords.push_back(MeshVertex[3 * vn]);
							pCoords.push_back(MeshVertex[3 * vn + 1]);
							pCoords.push_back(MeshVertex[3 * vn + 2]);
						}
						// add the new index to part's array
						pTris.push_back(pVert[vn]);
						// add the new index to temporary facet (for main array)
						vr[v] = mVert[vn];
					}
				}

				// save the temporary facet into main array in reverse order of vertices
				for (int v = 2; v >=0; v--) {
					mTris.push_back(vr[v]);
				}
			}
			
			// all the vertices is of the main
			else {
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in array yet
					if (mVert[vn] == -1) {
						// add the relation between solid's and main indices
						mVert[vn] = mvert++;
						// add the coordinates to main array
						mCoords.push_back(MeshVertex[3 * vn]);
						mCoords.push_back(MeshVertex[3 * vn + 1]);
						mCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to main array
					mTris.push_back(mVert[vn]);
				}
			}
		}

		try {
			// create the new CGAL mesh objects of new part and main
			Polyhedron part = arrToMesh(pCoords, pTris);
			Polyhedron main = arrToMesh(mCoords, mTris);

			// exclude the previous main mesh from widget
			tmesh.pop_back();
			// add the mesh of new part
			tmesh.push_back(part);
			// add the mesh of new main
			tmesh.push_back(main);
		}
		
		// the extracted vertices and facets don't form a proper mesh
		catch (const CGAL::Assertion_exception e) {
			// create the error message
			QMessageBox msg;
			msg.setText("The part extraction was unsuccessful, please reopen the model and try again");
			msg.setStandardButtons(QMessageBox::Ok);
			msg.setDefaultButton(QMessageBox::Ok);
			msg.exec();
		}
		// set off the previous skeleton
		vskel = 0;
		// recalculate arrays for OpenGL
		getMesh();
		getParts();

		return;
	}

	// we have no 'outer' branches; try to extract the branch from 'multi-thread' construction
	std::vector<uint32_t> seg3 = GetSeg3(vskel, eskel, SkelSegEdge, SkelSegVert);
	// check do we find such one
	if (seg3.size() > 0) {

		// initiate the arrays of new extracted part and of the rest of solid
		std::vector<int> pVert;			// part's vertices (relation to old indices)
		std::vector<double> pCoords;	// part's vertices coordinates
		std::vector<int> pTris;			// part's facets
		std::vector<int> mVert;			// main vertices (relation to old indices)
		std::vector<double> mCoords;	// main vertices coordinates
		std::vector<int> mTris;			// main facets

		// loop over solid's vertices + 2 new vertices from skeleton
		for (int j = 0; j < vmesh + 2; j++) {
			// initiate the arrays of indices relations
			pVert.push_back(-1);
			mVert.push_back(-1);
		}

		// create the new vertices from out vertices of skeleton's branch which will be used to cut solid
		// the start
		for (int v = 0; v < 3; v++) {
			// add it to both part's and main arrays
			pCoords.push_back(SkelVertex[3 * seg3[1] + v]);
			mCoords.push_back(SkelVertex[3 * seg3[1] + v]);
		}
		// the end
		for (int v = 0; v < 3; v++) {
			// add it to both part's and main arrays
			pCoords.push_back(SkelVertex[3 * seg3[2] + v]);
			mCoords.push_back(SkelVertex[3 * seg3[2] + v]);
		}
		// initiate the counters of number of vertices (part's and main)
		int pvert = 2;
		int mvert = 2;
		
		// loop over solid's facets
		for (int j = 0; j < fmesh; j++)
		{
			// initiate the checkers
			int isSeg = 0;		// number of vertices which is of new part
			int notSeg = -1;	// index of vertex which is not of the new part
			int isEnd = 0;		// the start/end checker

			// loop over facet's vertices
			for (int v = 0; v < 3; v++) {
				// set the checkers
				int skl = SkelMapMesh[MeshIndex[3 * j + v]];
				if (SkelSegVert[skl] == seg3[0] && skl != seg3[1] && skl != seg3[2]) isSeg++;
				else notSeg = v;
				if (skl == seg3[2] || skl == seg3[3]) isEnd = 1;
			}

			// all the vertices is of the new part
			if (isSeg == 3) {
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in array yet
					if (pVert[vn] == -1)
					{
						// add the relation between solid's and part's indices
						pVert[vn] = pvert++;
						// add the coordinates to part's array
						pCoords.push_back(MeshVertex[3 * vn]);
						pCoords.push_back(MeshVertex[3 * vn + 1]);
						pCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to part's array
					pTris.push_back(pVert[vn]);
				}
			}

			// the facet is on the border between part and main
			else if (isSeg == 2)
			{
				// save the existing faset to main
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in main array yet
					if (mVert[vn] == -1) {
						// add the relation between solid's and main indices
						mVert[vn] = mvert++;
						// add the coordinates to main array
						mCoords.push_back(MeshVertex[3 * vn]);
						mCoords.push_back(MeshVertex[3 * vn + 1]);
						mCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to main array
					mTris.push_back(mVert[vn]);
				}
				
				// create new facets to main and part
				// initiate the temporary facet
				std::array<int, 3> vr;
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					
					// change the vertex which is not of new part to skeleton's vertex
					if (notSeg == v) {
						// into part's array
						pTris.push_back(isEnd);
						// into temporary facet (for main array)
						vr[v] = isEnd;
					}
					// the vertices which is on new part will stay
					else {
						// it was not saved in part's array yet
						if (pVert[vn] == -1) {
							// add the relation between solid's and part's indices
							pVert[vn] = pvert++;
							// add the coordinates to part's array
							pCoords.push_back(MeshVertex[3 * vn]);
							pCoords.push_back(MeshVertex[3 * vn + 1]);
							pCoords.push_back(MeshVertex[3 * vn + 2]);
						}
						// add the new index to part's array
						pTris.push_back(pVert[vn]);
						// add the new index to temporary facet (for main array)
						vr[v] = mVert[vn];
					}
				}

				// save the temporary facet into main array in reverse order of vertices
				for (int v = 2; v >= 0; v--) {
					mTris.push_back(vr[v]);
				}
			}

			// all the vertices is of the main
			else {
				// loop over facet's vertices
				for (int v = 0; v < 3; v++) {
					// get the index of current vertex
					int vn = MeshIndex[3 * j + v];
					// it was not saved in array yet
					if (mVert[vn] == -1) {
						// add the relation between solid's and main indices
						mVert[vn] = mvert++;
						// add the coordinates to main array
						mCoords.push_back(MeshVertex[3 * vn]);
						mCoords.push_back(MeshVertex[3 * vn + 1]);
						mCoords.push_back(MeshVertex[3 * vn + 2]);
					}
					// add the new index to main array
					mTris.push_back(mVert[vn]);
				}
			}
		}

		try {
			// create the new CGAL mesh objects of new part and main
			Polyhedron part = arrToMesh(pCoords, pTris);
			Polyhedron main = arrToMesh(mCoords, mTris);

			// exclude the previous main mesh from widget
			tmesh.pop_back();
			// add the mesh of new part
			tmesh.push_back(part);
			// add the mesh of new main
			tmesh.push_back(main);
		}
		
		// the extracted vertices and facets don't form a proper mesh
		catch (const CGAL::Assertion_exception e) {
			// create the error message
			QMessageBox msg;
			msg.setText("The part extraction was unsuccessful, please reopen the model and try again");
			msg.setStandardButtons(QMessageBox::Ok);
			msg.setDefaultButton(QMessageBox::Ok);
			msg.exec();
		}

		// set off the previous skeleton
		vskel = 0;
		// recalculate arrays for OpenGL
		getMesh();
		getParts();

		return;
	}
	
	// we have no 'output' and 'multi-thread' branches, create the info message
	QMessageBox msg;
	msg.setText("Any acceptable type of part for extraction not found");
	msg.setStandardButtons(QMessageBox::Ok);
	msg.setDefaultButton(QMessageBox::Ok);
	msg.exec();

}

// Save the group of parts as separate OFF files
void Scene3D::groupsToOff()
{
	// get the filename using the standard Qt file dialog
	QString partName = QFileDialog::getSaveFileName(this, "Save", QDir::currentPath(), "Part (*.off)");
	// loop over the parts
	for (int part = 0; part < maxgroup + 1; part++) {
		
		// create the name of current part including the part index
		std::stringstream sc;
		sc.str(std::string());
		sc << std::setw(3) << std::setfill('0') << part << ".off";
		QString partN = QString::fromStdString(sc.str());

		// initiate the arrays 
		std::vector<int> NewVert;		// new indices of vertices
		std::vector<double> coords;		// vertices coodinates
		std::vector<int> tris;			// facets
		
		// loop over solid's vertices
		for (int i = 0; i < vmesh; i++) {
			
			// if the vertex is of target part
			if (SegGroups[SkelSegVert[SkelMapMesh[i]]] == part) {
				// add its coordinates to array
				coords.push_back(MeshVertex[3 * i]);
				coords.push_back(MeshVertex[3 * i + 1]);
				coords.push_back(MeshVertex[3 * i + 2]);
			}
			// add its new index to array
			NewVert.push_back((int)(coords.size() / 3) - 1);
		}

		// loop over facets
		for (int i = 0; i < fmesh; i++)
			// if the facet is of target part
			if (SegGroups[SkelSegVert[SkelMapMesh[MeshIndex[3 * i]]]] == part)
			{
				// add it to array
				tris.push_back(NewVert[MeshIndex[3 * i]]);
				tris.push_back(NewVert[MeshIndex[3 * i + 1]]);
				tris.push_back(NewVert[MeshIndex[3 * i + 2]]);
			}
		
		// convert the arrays into new CGAL mesh object
		Polyhedron partMesh = arrToMesh(coords, tris);
		// create the whole filename of current part
		QByteArray ba = partName.toLatin1();
		QString curName = QString(ba).replace(".off", partN);
		// save the part as OFF
		meshToOff(curName.toLatin1().data(), partMesh);
	}
}

// Draw the skeleton using Qt
void Scene3D::drawSkeleton()
{
	// check do the skeleton exist
	if (vskel == 0) return;
	// set the line width
	glLineWidth(5.0f);
	// set the vertices
	glVertexPointer(3, GL_FLOAT, 0, SkelVertex);
	// set the colors
	glColorPointer(4, GL_FLOAT, 0, SkelColor);
	// set the edges
	glDrawElements(GL_LINES, 2 * eskel, GL_UNSIGNED_INT, SkelIndex);
}

// Draw the facets of mesh
void Scene3D::drawFacets()
{
	if (mParentWnd->getSelectedSegment() == NULL) {
		// check do the mesh exist
		if (vmesh == 0) return;
		// set the vertices
		glVertexPointer(3, GL_FLOAT, 0, MeshVertex);
		// set the colors
		glColorPointer(4, GL_FLOAT, 0, MeshColor);
		// set the facets
		glDrawElements(GL_TRIANGLES, 3 * fmesh, GL_UNSIGNED_INT, MeshIndex);
	}
}

// Draw the wireframe of mesh
void Scene3D::drawWireframe()
{
	// check do the mesh exist
	if (vmesh == 0) return;
	// set the line width
	glLineWidth(2.0f);
	// set the vertices
	glVertexPointer(3, GL_FLOAT, 0, MeshVertex);
	// set the colors
	glColorPointer(4, GL_FLOAT, 0, EdgeColor);
	// set the edges
	glDrawElements(GL_LINES, 6 * fmesh, GL_UNSIGNED_INT, EdgeIndex);
}

// Draw the facets of extracted segments
void Scene3D::drawParts()
{
	// check do the segments exist
	if (vpart == 0) return;
	// set the vertices
	glVertexPointer(3, GL_FLOAT, 0, PartVertex);
	// set the colors
	glColorPointer(4, GL_FLOAT, 0, PartColor);
	// set the facets
	glDrawElements(GL_TRIANGLES, 3 * fpart, GL_UNSIGNED_INT, PartIndex);
}

// Draw the edges of extracted segments
void Scene3D::drawPWireframe()
{
	// check do the segments exist
	if (vpart == 0) return;
	// set the line width
	glLineWidth(2.0f);
	// set the vertices
	glVertexPointer(3, GL_FLOAT, 0, PartVertex);
	// set the colors
	glColorPointer(4, GL_FLOAT, 0, PEdgeColor);
	// set the edges
	glDrawElements(GL_LINES, 6 * fpart, GL_UNSIGNED_INT, PEdgeIndex);
}

// Initiate the scene with a new mesh
void Scene3D::load(Polyhedron mesh) 
{
	segMesh.clear();
	mParentWnd->clearTree();


	// loop over existing parts
	for (int i = (int)tmesh.size() - 1; i >= 0; i--) {
		// delete the current part
		Polyhedron p = tmesh[i];
		tmesh.pop_back();
		delete &p;
	}
	// add the input mesh into the scene's list
	tmesh.push_back(mesh);
	// reset the key index
	vskel = 0;
	// calculate the aspect ratio
	getAsp();
	// calculate the mesh arrays
	getMesh();
}


void Scene3D::updateMesh() {
	getParts();
	updateGL();
}


// Add new part to existing scene
void Scene3D::add(Polyhedron mesh)
{
	// add the input mesh into the scene's list
	tmesh.push_back(mesh);
	// reset the key index
	vskel = 0;
	// calculate the aspect ratio
	getAsp();
	// calculate the mesh arrays
	getMesh();
	// calculate the parts arrays
	getParts();
}


int Scene3D::testSegmentationBySDF() {
	int ret = EXIT_FAILURE;
	if (tmesh.size() == 1) {
		Polyhedron mesh = tmesh[0];

		// Segmentation via sdf values
		//// create a property-map for segment-ids
		//typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
		//Facet_int_map internal_segment_map;
		//boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

		//// calculate SDF values and segment the mesh using default parameters.
		//std::size_t number_of_segments = CGAL::segmentation_via_sdf_values(mesh, segment_property_map);

		//std::cout << "Number of segments: " << number_of_segments << std::endl;



		// Segmentation from sdf values
		// create a property-map for SDF values
		typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
		Facet_double_map internal_sdf_map;
		boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);
		// compute SDF values using default parameters for number of rays, and cone angle


		CGAL::sdf_values(	mesh, 
							sdf_property_map,
							2.0 / 3.0 * CGAL_PI,	// double cone_angle = 2.0 / 3.0 * CGAL_PI
							25,						// std::size_t number_of_rays = 25
							true					// bool postprocess = true
						);


		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
			//std::cout << sdf_property_map[facet_it] << " ";
			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			double value = sdf_property_map[facet_it];
			if (0 <= value && value <= 1.0 / 8.0) {
				R = 0.0;
				G = 0.0;
				B = 4.0 * value + .5; // .5 - 1 // b = 1/2
			}
			else if (1.0 / 8.0 < value && value <= 3.0 / 8.0) {
				R = 0.0;
				G = 4.0 * value - .5; // 0 - 1 // b = - 1/2
				B = 1.0; // small fix
			}
			else if (3.0 / 8.0 < value && value <= 5.0 / 8.0) {
				R = 4.0 * value - 1.5; // 0 - 1 // b = - 3/2
				G = 1.0;
				B = -4.0 * value + 2.5; // 1 - 0 // b = 5/2
			}
			else if (5.0 / 8.0 < value && value <= 7.0 / 8.0) {
				R = 1.0;
				G = -4.0 * value + 3.5; // 1 - 0 // b = 7/2
				B = 0.0;
			}
			else if (7.0 / 8.0 < value && value <= 1.0) {
				R = -4.0 * value + 4.5; // 1 - .5 // b = 9/2
				G = 0.0;
				B = 0.0;
			}
			else {    // should never happen - value > 1
				R = .5;
				G = 0;
				B = 0;
			}
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSDF[4 * vIndex] = R;
				MeshColorSDF[4 * vIndex + 1] = G;
				MeshColorSDF[4 * vIndex + 2] = B;
				MeshColorSDF[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());
		}



		// create a property-map for segment-ids
		typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
		Facet_int_map internal_segment_map;
		boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);

		// segment the mesh using default parameters for number of levels, and smoothing lambda
		// Any other scalar values can be used instead of using SDF values computed using the CGAL function


		std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(	mesh,
																				sdf_property_map,
																				segment_property_map,
																				5,							// std::size_t number_of_clusters = 5
																				0.26,						// double smoothing_lambda = 0.26
																				false						// bool output_cluster_ids = false
																			);
		std::cout << "Number of segments: " << number_of_segments << std::endl;




		std::vector<GLfloat> segmColorsArray;// [3 * number_of_segments];
		segmColorsArray.resize(3 * number_of_segments);
		// loop over segments
		for (int i = 0; i < number_of_segments; i++) {
			int mainCol = i % 3;
			// fill the colors with random (except of blue)
			segmColorsArray[3 * i] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 1] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + 2] = 0.05f + 0.01f * (qrand() % 15);
			segmColorsArray[3 * i + mainCol] += 0.75f;
		}


		// print segment-ids

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
			facet_it != mesh.facets_end(); ++facet_it) {
			std::cout << segment_property_map[facet_it] << " ";

			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			// set the different colors to segments if checked
			R = segmColorsArray[3 * segment_property_map[facet_it]];
			G = segmColorsArray[3 * segment_property_map[facet_it] + 1];
			B = segmColorsArray[3 * segment_property_map[facet_it] + 2];
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSeg[4 * vIndex] = R;
				MeshColorSeg[4 * vIndex + 1] = G;
				MeshColorSeg[4 * vIndex + 2] = B;
				MeshColorSeg[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());

		}
		std::cout << std::endl;

		ret = EXIT_SUCCESS;
	}

	return ret;
}

int Scene3D::testSegmentationBySkeletonAndSDF() {
	int ret = EXIT_FAILURE;
	if (tmesh.size() == 1) {
		Polyhedron mesh = tmesh[0];
		
		// extract the skeleton
		Skeleton skeleton;
		CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
		// init the polyhedron simplex indices
		CGAL::set_halfedgeds_items_id(mesh);
		//for each input vertex compute its distance to the skeleton
		std::vector<double> distances(num_vertices(mesh));
		BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton))
		{
			const Point& skel_pt = skeleton[v].point;
			BOOST_FOREACH(vertex_descriptor mesh_v, skeleton[v].vertices)
			{
				const Point& mesh_pt = mesh_v->point();
				distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
			}
		}

		// create a property-map for sdf values
		std::vector<double> sdf_values(num_faces(mesh));
		Facet_with_id_pmap<double> sdf_property_map(sdf_values);
		// compute sdf values with skeleton
		BOOST_FOREACH(face_descriptor f, faces(mesh))
		{
			double dist = 0;
			BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(halfedge(f, mesh), mesh))
				dist += distances[target(hd, mesh)->id()];
			sdf_property_map[f] = dist / 3.;
		}
		// post-process the sdf values
		CGAL::sdf_values_postprocessing(mesh, sdf_property_map);

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
			//std::cout << sdf_property_map[facet_it] << " ";
			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			double value = sdf_property_map[facet_it];
			if (0 <= value && value <= 1.0 / 8.0) {
				R = 0.0;
				G = 0.0;
				B = 4.0 * value + .5; // .5 - 1 // b = 1/2
			}
			else if (1.0 / 8.0 < value && value <= 3.0 / 8.0) {
				R = 0.0;
				G = 4.0 * value - .5; // 0 - 1 // b = - 1/2
				B = 1.0; // small fix
			}
			else if (3.0 / 8.0 < value && value <= 5.0 / 8.0) {
				R = 4.0 * value - 1.5; // 0 - 1 // b = - 3/2
				G = 1.0;
				B = -4.0 * value + 2.5; // 1 - 0 // b = 5/2
			}
			else if (5.0 / 8.0 < value && value <= 7.0 / 8.0) {
				R = 1.0;
				G = -4.0 * value + 3.5; // 1 - 0 // b = 7/2
				B = 0.0;
			}
			else if (7.0 / 8.0 < value && value <= 1.0) {
				R = -4.0 * value + 4.5; // 1 - .5 // b = 9/2
				G = 0.0;
				B = 0.0;
			}
			else {    // should never happen - value > 1
				R = .5;
				G = 0;
				B = 0;
			}
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSDF[4 * vIndex] = R;
				MeshColorSDF[4 * vIndex + 1] = G;
				MeshColorSDF[4 * vIndex + 2] = B;
				MeshColorSDF[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());
		}

		// create a property-map for segment-ids (it is an adaptor for this case)
		std::vector<std::size_t> segment_ids(num_faces(mesh));
		Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
		// segment the mesh using default parameters
		std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map);


		if (SkelSegmColors != NULL) {
			delete SkelSegmColors;
		}
		SkelSegmColors = new GLfloat[3 * number_of_segments];
		// loop over segments
		for (int i = 0; i < number_of_segments; i++) {
			int mainCol = i % 3;
			// fill the colors with random (except of blue)
			SkelSegmColors[3 * i] = 0.05f + 0.01f * (qrand() % 15);
			SkelSegmColors[3 * i + 1] = 0.05f + 0.01f * (qrand() % 15);
			SkelSegmColors[3 * i + 2] = 0.05f + 0.01f * (qrand() % 15);
			SkelSegmColors[3 * i + mainCol] += 0.75f;
		}


		// print segment-ids

		for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
			facet_it != mesh.facets_end(); ++facet_it) {
			std::cout << segment_property_map[facet_it] << " ";

			GLfloat R = 0.5f, G = 0.7f, B = 0.5f, T = 0.3f;
			// set the different colors to segments if checked
			R = SkelSegmColors[3 * segment_property_map[facet_it]];
			G = SkelSegmColors[3 * segment_property_map[facet_it] + 1];
			B = SkelSegmColors[3 * segment_property_map[facet_it] + 2];
			T = 1.0;

			Halfedge_facet_circulator j = facet_it->facet_begin();
			// Facets in polyhedral surfaces are at least triangles.
			CGAL_assertion(CGAL::circulator_size(j) >= 3);
			std::cout << CGAL::circulator_size(j) << ' ';
			do {
				int vIndex = std::distance(mesh.vertices_begin(), j->vertex());
				//std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());

				MeshColorSeg[4 * vIndex] = R;
				MeshColorSeg[4 * vIndex + 1] = G;
				MeshColorSeg[4 * vIndex + 2] = B;
				MeshColorSeg[4 * vIndex + 3] = T;

			} while (++j != facet_it->facet_begin());
		}
		ret = EXIT_SUCCESS;
	}
	return ret;
}

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel1;
typedef CGAL::Polyhedron_3<Kernel1> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel1, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

int Scene3D::testPolyedraDecomposition() {
	int ret = EXIT_FAILURE;
	//if (true/*tmesh.size() == 1*/) {
	//	Polyhedron_3 mesh;// = tmesh[0];
	//	std::ifstream inStrm("F:\\Lavoro\\MElmasry\\Part1.off", std::ios::in);
	//	inStrm >> mesh;

	//	Nef_polyhedron_3 N(mesh);
	//	CGAL::convex_decomposition_3(N);
	//	std::vector<Polyhedron_3> convex_parts;
	//	Volume_const_iterator ci = ++N.volumes_begin();
	//	for (; ci != N.volumes_end(); ++ci) {
	//		if (ci->mark()) {
	//			Polyhedron_3 P;
	//			N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
	//			convex_parts.push_back(P);
	//		}
	//	}
	//	std::string fName;
	//	for (int i = 0 ; i < (int)convex_parts.size() ; i++) {
	//		std::ostringstream stringStream;
	//		stringStream << "C:\\temp\\convex_part_";
	//		stringStream << i;
	//		stringStream << ".off";
	//		fName = stringStream.str();
	//		std::ofstream outStrm(fName.c_str(), std::ios::out);
	//		outStrm << convex_parts[i];
	//		outStrm.close();
	//	}
	//	ret = EXIT_SUCCESS;
	//}
	return ret;
}

int Scene3D::updateSegmentationTree() {
	int ret = EXIT_FAILURE;
	mParentWnd->clearTree();
	if (SegGraph.nodeCount() > 0) {
		ret = EXIT_SUCCESS;
		int segmentTextCounter = 1;
		for (int i = 0 ; i < SegGraph.hRootCount() ; i++) {

			mParentWnd->addTreeItem(SegGraph.getHRootAt(i), NULL, QString("Segment%1").arg(segmentTextCounter));
			segmentTextCounter++;

			std::stack<SegmentGraph::SegmentNode*> stackNode;
			stackNode.push(SegGraph.getHRootAt(i));
			while (stackNode.empty() == false) {
				SegmentGraph::SegmentNode *n = stackNode.top();
				stackNode.pop();
				for (int j = 0 ; j < n->childCount() ; j++) {
					SegmentGraph::SegmentNode *childNode = n->getChildAt(j);
					mParentWnd->addTreeItem(childNode, n, QString("Segment%1").arg(segmentTextCounter));
					segmentTextCounter++;
					stackNode.push(childNode);
				}
			}
		}
	}

	return ret;
}

int Scene3D::computeSegmentHierarchy() {
	int ret = EXIT_FAILURE;
	if (segMesh.empty() == false) {
		if (mParentWnd != NULL) {
			mParentWnd->clearTree();
		}
		SegGraph.clear();
		std::vector<CGAL::Bbox_3> segBBoxArray;
		for (int i = 0 ; i < (int)segMesh.size() ; i++) {
			SegGraph.addSegment(&segMesh[i]);
			CGAL::Bbox_3 bbox;
			segBBoxArray.push_back(bbox);
			for (Polyhedron::Facet_iterator facet_it = segMesh[i].facets_begin(); facet_it != segMesh[i].facets_end(); ++facet_it) {
				Halfedge_facet_circulator fj = facet_it->facet_begin();
				int j = 0;
				do {
					segBBoxArray[i] += fj->vertex()->point().bbox();
				} while (++fj != facet_it->facet_begin());
			}
		}

		for (int i = 0 ; i < SegGraph.nodeCount() ; i++) {
			for (int j = i + 1; j < SegGraph.nodeCount(); j++) {
				Polyhedron *poly1 = SegGraph.getNodeAt(i)->segment();
				Polyhedron *poly2 = SegGraph.getNodeAt(j)->segment();


				if (poly1 != NULL && poly2 != NULL && poly1 != poly2) {
					if (CGAL::Polygon_mesh_processing::do_intersect(*poly1, *poly2) == true) {
						SegGraph.addConnection(poly1, poly2);
					}
				}
			}
		}

		int minNodeIndex = -1;
		do {
			minNodeIndex = -1;
			double minZ = DBL_MAX;
			for (int i = 0 ; i < SegGraph.nodeCount() ; i++) {
				if (SegGraph.getNodeAt(i)->getFlag() == SEGGRAPH_NOTMARKED) {
					if (segBBoxArray[i].zmin() < minZ) {
						minZ = segBBoxArray[i].zmin();
						minNodeIndex = i;
					}
				}
			}
			if (minNodeIndex != -1) {
				SegmentGraph::SegmentNode *n = SegGraph.getNodeAt(minNodeIndex);
				if (n != NULL) {
					n->setFlag(SEGGRAPH_VISITED);
					SegGraph.addHRoot(n);

					std::stack<SegmentGraph::SegmentNode*> nodeStack;
					nodeStack.push(n);

					while (nodeStack.empty() == false) {
						SegmentGraph::SegmentNode *node = nodeStack.top();
						nodeStack.pop();
						std::vector<SegmentGraph::SegmentNode*> connectedNodes = SegGraph.getConnectedNodes(node);
						for (int i = 0 ; i < (int)connectedNodes.size() ; i++) {
							if (connectedNodes[i]->getFlag() == SEGGRAPH_NOTMARKED) {
								connectedNodes[i]->setFlag(SEGGRAPH_VISITED);
								node->appendChild(connectedNodes[i]);
								nodeStack.push(connectedNodes[i]);
							}
						}
					}
				}
			}
		} while (minNodeIndex != -1);

		ret = EXIT_SUCCESS;
	}
	return ret;
}





/* ---------------------------------------------------------------------------------------------------------- */


SegmentGraph::SegmentNode::SegmentNode() {
	mFlag = 0;
	pSegment = NULL;
	mChildList.clear();
}

SegmentGraph::SegmentNode::SegmentNode(Polyhedron *_seg) {
	mFlag = 0;
	pSegment = _seg;
	mChildList.clear();
}


SegmentGraph::SegmentNode::~SegmentNode() {
	mChildList.clear();
}


Polyhedron* SegmentGraph::SegmentNode::segment() const {
	return pSegment;
}

void SegmentGraph::SegmentNode::segment(Polyhedron* _seg) {
	pSegment = _seg;
}

int SegmentGraph::SegmentNode::childCount() const {
	return (int)mChildList.size();
}

SegmentGraph::SegmentNode* SegmentGraph::SegmentNode::getChildAt(int _index) const {
	SegmentNode *ret = NULL;
	if (_index >= 0 && _index < mChildList.size()) {
		ret = mChildList[_index];
	}
	return ret;
}

void SegmentGraph::SegmentNode::appendChild(SegmentGraph::SegmentNode* _node) {
	mChildList.push_back(_node);
}


void SegmentGraph::SegmentNode::insertChildAfter(SegmentGraph::SegmentNode* _node, int _index) {
	if (_index >= 0 && _index < mChildList.size()) {
		if (_index + 1 < mChildList.size()) {
			mChildList.insert(mChildList.begin() + _index + 1, _node);
		}
		else {
			mChildList.push_back(_node);
		}
	}
}



SegmentGraph::SegmentNode* SegmentGraph::SegmentNode::removeChildAt(int _index) {
	SegmentNode *ret = NULL;
	if (_index >= 0 && _index < mChildList.size()) {
		ret = mChildList[_index];
		mChildList.erase(mChildList.begin() + _index);
	}
	return ret;
}


void SegmentGraph::SegmentNode::clearChildList() {
	mChildList.clear();
}


int SegmentGraph::SegmentNode::getFlag() const {
	return mFlag;
}


void SegmentGraph::SegmentNode::setFlag(int _flag) {
	mFlag = _flag;
}


/* ---------------------------------------------------------------------------------------------------------- */



SegmentGraph::SegmentGraph() {
	mNodeList.clear();
	mEdgeList.clear();
	mHRootList.clear();
}


SegmentGraph::~SegmentGraph() {
	clear();
}



void SegmentGraph::addSegment(Polyhedron *_seg) {
	SegmentGraph::SegmentNode *node = new SegmentGraph::SegmentNode(_seg);
	mNodeList.push_back(node);
}


SegmentGraph::SegmentNode* SegmentGraph::getNodeWithSegment(Polyhedron *_seg) const {
	SegmentGraph::SegmentNode *ret = NULL;
	for (int i = 0; i < (int)mNodeList.size() && ret == NULL; i++) {
		if (mNodeList[i]->segment() == _seg) {
			ret = mNodeList[i];
		}
	}
	return ret;
}


void SegmentGraph::addConnection(Polyhedron *_seg1, Polyhedron *_seg2) {
	SegmentGraph::SegmentNode *n1 = getNodeWithSegment(_seg1);
	SegmentGraph::SegmentNode *n2 = getNodeWithSegment(_seg2);
	if (n1 != NULL && n2 != NULL) {
		mEdgeList.push_back(std::pair<SegmentGraph::SegmentNode*, SegmentGraph::SegmentNode*>(n1, n2));
	}
}

void SegmentGraph::addChild(Polyhedron *_father, Polyhedron *_child) {
	SegmentGraph::SegmentNode *fatherNode = getNodeWithSegment(_father);
	SegmentGraph::SegmentNode *childNode = getNodeWithSegment(_child);
	if (fatherNode != NULL && childNode != NULL) {
		fatherNode->appendChild(childNode);
	}
}


void SegmentGraph::addChild(SegmentGraph::SegmentNode *_father, SegmentGraph::SegmentNode *_child) {
	if (_father != NULL && _child != NULL) {
		_father->appendChild(_child);
	}
}


void SegmentGraph::addHRoot(Polyhedron *_seg) {
	SegmentGraph::SegmentNode *n = getNodeWithSegment(_seg);
	if (n != NULL) {
		mHRootList.push_back(n);
	}
}

void SegmentGraph::addHRoot(SegmentGraph::SegmentNode *_node) {
	if (_node != NULL) {
		mHRootList.push_back(_node);
	}
}


void SegmentGraph::clear() {
	for (int i = 0; i < (int)mNodeList.size(); i++) {
		if (mNodeList[i] != NULL) {
			delete mNodeList[i];
		}
	}
	mNodeList.clear();
	mEdgeList.clear();
	mHRootList.clear();
}


int SegmentGraph::nodeCount() const {
	return (int)mNodeList.size();
}


SegmentGraph::SegmentNode* SegmentGraph::getNodeAt(int _index) const {
	SegmentNode *ret = NULL;
	if (_index >= 0 && _index < mNodeList.size()) {
		ret = mNodeList[_index];
	}
	return ret;
}


std::vector<SegmentGraph::SegmentNode*> SegmentGraph::getConnectedNodes(SegmentGraph::SegmentNode *_node) const {
	std::vector<SegmentGraph::SegmentNode*> ret;
	for (int i = 0 ; i < (int)mEdgeList.size() ; i++) {
		if (mEdgeList[i].first == _node) {
			ret.push_back(mEdgeList[i].second);
		}
		else if (mEdgeList[i].second == _node) {
			ret.push_back(mEdgeList[i].first);
		}
	}
	return ret;
}


SegmentGraph::SegmentNode* SegmentGraph::getHRootAt(int _index) const {
	SegmentNode *ret = NULL;
	if (_index >= 0 && _index < mHRootList.size()) {
		ret = mHRootList[_index];
	}
	return ret;
}


int SegmentGraph::hRootCount() const {
	return (int)mHRootList.size();
}
