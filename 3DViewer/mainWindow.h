#ifndef MAINWINDOW_H
#define MAINWINDOW_H

// The masks of 'elements visibility' variable 
#define shAxis 0x01
#define shSkeleton 0x02
#define shWireframe 0x04
#define shFacets 0x08
#define shParts 0x10
#define shColors 0x20
#define segColors 0x40
#define sdfColors 0x80

#include <QMainWindow>
#include <qtreewidget.h>
#include "scene3d.h"

class SegmentGraph::SegmentNode;

class MainWindow : public QMainWindow
{
	Q_OBJECT

private:


public:
	MainWindow();

	std::vector<std::pair<QTreeWidgetItem*, SegmentGraph::SegmentNode*> > mTreeItems;

	//QList<QTreeWidgetItem*> mItems;
	QHBoxLayout *mLayout;
	QTreeWidget *mTree;
	Scene3D * widget;		// Qt widget to show the 3D objects
	QMenu * menuActions;	// 'Process' menu
	QMenu * menuOptions;	// 'Elements' menu
	QMenu * menuTest;

	QTreeWidgetItem *mSelectedTreeItem;

	void clearTree();
	void addTreeItem(SegmentGraph::SegmentNode *_node, SegmentGraph::SegmentNode *_parentNode, QString _text);
	SegmentGraph::SegmentNode* getSelectedSegment();

private slots:

	void openModel();
	void addPart();
	void savePOFF();
	void saveSOFF();
	void Skeleton();
	void sGroups();
	void nPart();
	int createSkeleton() { return widget->getSkeleton(); }
	void sColors() { widget->switchColors(); }
	void saveGroups() { widget->groupsToOff(); }
	void setFigureOn();

	void setDockOptions();
	void setDockOptions1();
	void setDockOptions2();
	void setDockOptions3();

	void keyPressEvent(QKeyEvent* pe) { return widget->keyPressEvent(pe); }
	int TestSegmentationBySDF();
	int TestSegmentationBySkeletonAndSDF();
	int TestPolyedraDecomposition();

	void splitAndSaveBySkeleton();
	void splitAndSaveBySDF();
	void splitAndSaveBySkeletonAndSDF();

	void segmentItemClicked(QTreeWidgetItem * item, int column);

};
#endif